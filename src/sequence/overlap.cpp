//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <thread>
#include <ctime>
#include <queue>
#include <cmath>
#include <chrono>
#include <ctime>
#include <cstring>
#include <iomanip>
#include <numeric>

#include "overlap.h"
#include "alignment.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "../common/disjoint_set.h"
#include "../common/bfcontainer.h"


//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp,
								  bool forceLocal) const
{
	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	//filter overlaps that way to divergent in length.
	//theoretically, they should not pass sequence divergence filter,
	//but just in case
	static const float OVLP_DIV = 0.5;
	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	if (lengthDiff > OVLP_DIV * std::min(ovlp.curRange(), ovlp.extRange()))
	{
		return false;
	}

	//check if it's "almost trivial" match with intersecting sequence
	if (ovlp.curId == ovlp.extId)
	{
		int32_t intersect = std::min(ovlp.curEnd, ovlp.extEnd) - 
			   				std::max(ovlp.curBegin, ovlp.extBegin);
		if (intersect > ovlp.curRange() / 2) return false;
	}

	//check "strand skipping" PacBio pattern
	if (ovlp.curId == ovlp.extId.rc())
	{
		int32_t intersect = std::min(ovlp.curEnd, ovlp.extLen - ovlp.extBegin) - 
			   				std::max(ovlp.curBegin, ovlp.extLen - ovlp.extEnd);

		//if (intersect > -_maxJump) outSuggestChimeric = true;
		if (intersect > ovlp.curRange() / 2) return false;
	}

	if (!forceLocal && _checkOverhang && ovlp.lrOverhang() > _maxOverhang) return false;

	return true;
}


namespace
{
	struct KmerMatch
	{
		KmerMatch(int32_t cur = 0, int32_t ext = 0,
				  FastaRecord::Id extId = FastaRecord::ID_NONE): 
			curPos(cur), extPos(ext), extId(extId) {}
		int32_t curPos;
		int32_t extPos;
		FastaRecord::Id extId;
	};

	template <class T>
	void shrinkAndClear(std::vector<T>& vec, float rate)
	{
		if (vec.empty()) return;
		//if ((float)vec.capacity() / vec.size() < rate) return;
		//Logger::get().debug() << "Shrink: " << vec.capacity() << " " << vec.size();

		size_t newCapacity = std::roundf((float)vec.capacity() / rate);
		vec = std::vector<T>();
		vec.reserve(newCapacity);
	}
}

//This implementation was inspired by Heng Li's minimap2 paper
//might be used in parallel
std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec, 
								bool forceLocal,
								OvlpDivStats& divStats,
								int maxOverlaps) const
{
	//static std::ofstream fout("../kmers.txt");
	
	//const int MAX_LOOK_BACK = 50;
	const int kmerSize = Parameters::get().kmerSize;
	//const float minKmerSruvivalRate = std::exp(-_maxDivergence * kmerSize);
	const float minKmerSruvivalRate = 0.01;

	static const float LG_GAP = (float)Config::get("chain_large_gap_penalty");
	static const float SM_GAP = (float)Config::get("chain_small_gap_penalty");
	static const int GAP_JUMP_THLD = (int)Config::get("chain_gap_jump_threshold");

	//outSuggestChimeric = false;
	int32_t curLen = fastaRec.sequence.length();
	std::vector<int32_t> curFilteredPos;

	//cache memory-intensive containers as
	//many parallel memory allocations slow us down significantly
	//thread_local std::vector<KmerMatch> vecMatches;
	thread_local std::vector<KmerMatch> matchesList;
	thread_local std::vector<int32_t> scoreTable;
	thread_local std::vector<int32_t> backtrackTable;

	static ChunkPool<KmerMatch> sharedChunkPool;	//shared accoress threads
	BFContainer<KmerMatch> vecMatches(sharedChunkPool);

	//speed benchmarks
	thread_local float timeMemory = 0;
	thread_local float timeKmerIndexFirst = 0;
	thread_local float timeKmerIndexSecond = 0;
	thread_local float timeDp = 0;
	auto timeStart = std::chrono::system_clock::now();

	/*static std::mutex reportMutex;
	thread_local int threadId = rand() % 1000; 
	thread_local int numTicks = 0;
	++numTicks;
	thread_local auto prevReport = std::chrono::system_clock::now();
	if ((std::chrono::system_clock::now() - prevReport) > 
		std::chrono::seconds(60))
	{
		std::lock_guard<std::mutex> lock(reportMutex);
		prevReport = std::chrono::system_clock::now();
		Logger::get().debug() << ">Perf " << threadId
			<< " mem:" << timeMemory
			<< " kmFst:" << timeKmerIndexFirst 
			<< " kmSnd:" << timeKmerIndexSecond 
			<< " dp:" << timeDp << " ticks: " << numTicks; 
		Logger::get().debug() << ">Mem  " << threadId 
			<< " chunks:" << sharedChunkPool.numberChunks();

		timeMemory = 0;
		timeKmerIndexFirst = 0;
		timeKmerIndexSecond = 0;
		timeDp = 0;
		numTicks = 0;
	}*/
	timeStart = std::chrono::system_clock::now();

	//although once in a while shrink allocated memory size
	//thread_local auto prevCleanup = 
	//	std::chrono::system_clock::now() + std::chrono::seconds(rand() % 60);
	thread_local int prevCleanup = 0;
	if (++prevCleanup > 50)
	{
		prevCleanup = 0;
		shrinkAndClear(matchesList, 2);
		shrinkAndClear(scoreTable, 2);
		shrinkAndClear(backtrackTable, 2);
	}
	timeMemory += std::chrono::duration_cast<std::chrono::duration<float>>
						(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	for (const auto& curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (_vertexIndex.isRepetitive(curKmerPos.kmer))
		{
			curFilteredPos.push_back(curKmerPos.position);
			continue;
		}
		if (!_vertexIndex.kmerFreq(curKmerPos.kmer)) continue;
		if (!_vertexIndex.isWhitelisted(curKmerPos.kmer)) continue;

		//FastaRecord::Id prevSeqId = FastaRecord::ID_NONE;
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(curKmerPos.kmer))
		{
			//no trivial matches
			if ((extReadPos.readId == fastaRec.id &&
				extReadPos.position == curKmerPos.position)) continue;

			vecMatches.emplace_back(curKmerPos.position, 
									extReadPos.position,
									extReadPos.readId);
		}
	}
	timeKmerIndexFirst += std::chrono::duration_cast<std::chrono::duration<float>>
							(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	std::sort(vecMatches.begin(), vecMatches.end(),
			  [](const KmerMatch& k1, const KmerMatch& k2)
			  {return k1.extId != k2.extId ? k1.extId < k2.extId : 
			  								 k1.curPos < k2.curPos;});

	timeKmerIndexSecond += std::chrono::duration_cast<std::chrono::duration<float>>
								(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	const int STAT_WND = 10000;
	std::vector<OverlapRange> divStatWindows(curLen / STAT_WND + 1);

	std::vector<OverlapRange> detectedOverlaps;
	size_t extRangeBegin = 0;
	size_t extRangeEnd = 0;
	while(extRangeEnd < vecMatches.size())
	{
		if (maxOverlaps != 0 &&
			detectedOverlaps.size() >= (size_t)maxOverlaps) break;

		extRangeBegin = extRangeEnd;
		size_t uniqueMatches = 0;
		int32_t prevPos = 0;
		while (extRangeEnd < vecMatches.size() &&
			   vecMatches[extRangeBegin].extId == 
			   vecMatches[extRangeEnd].extId)
		{
			if (vecMatches[extRangeEnd].curPos != prevPos)
			{
				++uniqueMatches;
				prevPos = vecMatches[extRangeEnd].curPos;
			}
			++extRangeEnd;
		}
		if (uniqueMatches < minKmerSruvivalRate * _minOverlap) continue;

		matchesList.assign(vecMatches.begin() + extRangeBegin,
						   vecMatches.begin() + extRangeEnd);
		assert(matchesList.size() > 0 && 
			   matchesList.size() < (size_t)std::numeric_limits<int32_t>::max());

		FastaRecord::Id extId = matchesList.front().extId;
		int32_t extLen = _seqContainer.seqLen(extId);

		//pre-filtering
		int32_t minCur = matchesList.front().curPos;
		int32_t maxCur = matchesList.back().curPos;
		int32_t minExt = std::numeric_limits<int32_t>::max();
		int32_t maxExt = std::numeric_limits<int32_t>::min();
		for (const auto& match : matchesList)
		{
			minExt = std::min(minExt, match.extPos);
			maxExt = std::max(maxExt, match.extPos);
		}
		if (maxCur - minCur < _minOverlap || 
			maxExt - minExt < _minOverlap) continue;
		if (_checkOverhang && !forceLocal)
		{
			if (std::min(minCur, minExt) > _maxOverhang) continue;
			if (std::min(curLen - maxCur, 
						 extLen - maxExt) > _maxOverhang) continue;
		}
		//++uniqueCandidates;

		//chain matiching positions with DP
		scoreTable.assign(matchesList.size(), 0);
		backtrackTable.assign(matchesList.size(), -1);

		bool extSorted = extLen > curLen;
		if (extSorted)
		{
			std::sort(matchesList.begin(), matchesList.end(),
					  [](const KmerMatch& k1, const KmerMatch& k2)
					  {return k1.extPos < k2.extPos;});
		}

		for (int32_t i = 1; i < (int32_t)scoreTable.size(); ++i)
		{
			int32_t maxScore = 0;
			int32_t maxId = 0;
			int32_t curNext = matchesList[i].curPos;
			int32_t extNext = matchesList[i].extPos;
			//int32_t noImprovement = 0;

			for (int32_t j = i - 1; j >= 0; --j)
			{
				int32_t curPrev = matchesList[j].curPos;
				int32_t extPrev = matchesList[j].extPos;
				if (0 < curNext - curPrev && curNext - curPrev < _maxJump &&
					0 < extNext - extPrev && extNext - extPrev < _maxJump)
				{
					int32_t matchScore = 
						std::min(std::min(curNext - curPrev, extNext - extPrev),
										  kmerSize);
					int32_t jumpDiv = abs((curNext - curPrev) - 
										  (extNext - extPrev));
					//int32_t gapCost = jumpDiv ? 
					//		kmerSize * jumpDiv + ilog2_32(jumpDiv) : 0;
					int32_t gapCost = (jumpDiv > GAP_JUMP_THLD ? LG_GAP : SM_GAP) * jumpDiv;
					int32_t nextScore = scoreTable[j] + matchScore - gapCost;
					if (nextScore > maxScore)
					{
						maxScore = nextScore;
						maxId = j;
						//noImprovement = 0;

						if (jumpDiv == 0 && curNext - curPrev < kmerSize) break;
					}
					/*else
					{
						if (++noImprovement > MAX_LOOK_BACK) break;
					}*/
				}
				if (extSorted && extNext - extPrev > _maxJump) break;
				if (!extSorted && curNext - curPrev > _maxJump) break;
			}

			scoreTable[i] = std::max(maxScore, kmerSize);
			if (maxScore > kmerSize)
			{
				backtrackTable[i] = maxId;
			}
		}

		//backtracking
		std::vector<OverlapRange> extOverlaps;
		//std::vector<int32_t> shifts;
		std::vector<std::pair<int32_t, int32_t>> kmerMatches;
		
		//initiate chains from the highest scores in the table (e.g. local maximums)
		std::vector<size_t> orderedScores(backtrackTable.size());
		std::iota(orderedScores.begin(), orderedScores.end(), 0);
		std::sort(orderedScores.begin(), orderedScores.end(),
				  [](size_t a, size_t b) {return scoreTable[a] > scoreTable[b];});

		//for (int32_t chainStart = backtrackTable.size() - 1; 
		//	 chainStart > 0; --chainStart)
		for (int32_t chainStart : orderedScores)
		{
			if (backtrackTable[chainStart] == -1) continue;

			//int32_t chainMaxScore = scoreTable[chainStart];
			int32_t lastMatch = chainStart;
			int32_t firstMatch = 0;

			int32_t chainLength = 0;
			//shifts.clear();
			kmerMatches.clear();
			
			int32_t pos = chainStart;
			while (pos != -1)
			{
				//found a new maximum, shorten the chain end
				/*if (scoreTable[pos] > chainMaxScore)
				{
					chainMaxScore = scoreTable[pos];
					lastMatch = pos;
					chainLength = 0;
					shifts.clear();
					kmerMatches.clear();
				}*/

				firstMatch = pos;
				//shifts.push_back(matchesList[pos].curPos - 
				//				 matchesList[pos].extPos);
				++chainLength;

				if (_keepAlignment)
				{
					if (kmerMatches.empty() || 
						kmerMatches.back().first - matchesList[pos].curPos >
						kmerSize)
					{
						kmerMatches.emplace_back(matchesList[pos].curPos,
								 				 matchesList[pos].extPos);
					}
				}

				assert(pos >= 0 && pos < (int32_t)backtrackTable.size());
				int32_t newPos = backtrackTable[pos];
				backtrackTable[pos] = -1;
				pos = newPos;
			}

			//Logger::get().debug() << chainStart - firstMatch << " " << lastMatch - firstMatch;

			OverlapRange ovlp(fastaRec.id, matchesList.front().extId,
							  matchesList[firstMatch].curPos, 
							  matchesList[firstMatch].extPos,
							  curLen, extLen);
			ovlp.curEnd = matchesList[lastMatch].curPos + kmerSize - 1;
			ovlp.extEnd = matchesList[lastMatch].extPos + kmerSize - 1;
			ovlp.score = scoreTable[lastMatch] - scoreTable[firstMatch] + 
						 kmerSize - 1;

			if (this->overlapTest(ovlp, forceLocal))
			{
				if (_keepAlignment)
				{
					kmerMatches.emplace_back(ovlp.curBegin, ovlp.extBegin);
					std::reverse(kmerMatches.begin(), kmerMatches.end());
					kmerMatches.emplace_back(ovlp.curEnd, ovlp.extEnd);
					ovlp.kmerMatches = new std::vector<std::pair<int32_t, int32_t>>();
					ovlp.kmerMatches->swap(kmerMatches);
				}
				//ovlp.leftShift = median(shifts);
				//ovlp.rightShift = extLen - curLen + ovlp.leftShift;

				//estimating identity using k-mers
				int32_t filteredPositions = 0;
				for (auto pos : curFilteredPos)
				{
					if (pos < ovlp.curBegin) continue;
					if (pos > ovlp.curEnd) break;
					++filteredPositions;
				}
				float normLen = std::max(ovlp.curRange(), 
										 ovlp.extRange()) - filteredPositions;
				float matchRate = (float)chainLength * 
								  _vertexIndex.getSampleRate() / normLen;
				matchRate = std::min(matchRate, 1.0f);
				//float repeatRate = (float)filteredPositions / ovlp.curRange();
				ovlp.seqDivergence = std::log(1 / matchRate) / kmerSize;
				//ovlp.seqDivergence += _estimatorBias;
				extOverlaps.push_back(ovlp);
			}
		}

		//now we have a lits of (possibly multiple) putative overlaps
		//agains a singe ext sequence. Now, select the list of primary overlaps
		std::vector<OverlapRange> primaryOverlaps;
		std::sort(extOverlaps.begin(), extOverlaps.end(),
				  [](const OverlapRange& r1, const OverlapRange& r2)
				  {return r1.score > r2.score;});

		if (_onlyMaxExt)	//select single best overlap
		{
			if (!extOverlaps.empty()) primaryOverlaps.push_back(extOverlaps.front());
		}
		else
		{
			for (const auto& ovlp : extOverlaps)
			{
				bool isContained = false;
				for (const auto& prim : primaryOverlaps)
				{
					if (ovlp.containedBy(prim) && prim.score > ovlp.score)
					{
						isContained = true;
						break;
					}
				}
				if (!isContained)
				{
					primaryOverlaps.push_back(ovlp);
				}
			}
		}

		//divergence check for the selected primary overlaps
		for (auto& ovlp : primaryOverlaps)
		{
			if(_nuclAlignment)	//identity using base-level alignment
			{
				ovlp.seqDivergence = getAlignmentErrEdlib(ovlp, fastaRec.sequence, 
														   _seqContainer.getSeq(extId),
														   _maxDivergence, _useHpc);
			}

			if (ovlp.seqDivergence < _maxDivergence)
			{
				detectedOverlaps.push_back(ovlp);
			}
			//if alignment not passing thrshold, check if its parts do
			else if (_partitionBadMappings)
			{
				auto trimmedOverlaps = 
					checkIdyAndTrim(ovlp, fastaRec.sequence, 
								    _seqContainer.getSeq(extId),
								    _maxDivergence, _minOverlap, _useHpc);
				for (auto& trimOvlp : trimmedOverlaps)
				{
					detectedOverlaps.push_back(trimOvlp);
				}
			}

			//statistics
			size_t wnd = ovlp.curBegin / STAT_WND;
			if (ovlp.curRange() > divStatWindows[wnd].curRange())
			{
				divStatWindows[wnd] = ovlp;
			}
		}
	}

	timeDp += std::chrono::duration_cast<std::chrono::duration<float>>
				(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	for (const auto& ovlp : divStatWindows)
	{
		if (ovlp.curRange() > 0)
		{
			divStats.add(ovlp.seqDivergence);
		}
	}
	return detectedOverlaps;
}

bool OverlapContainer::hasSelfOverlaps(FastaRecord::Id readId)
{
	this->lazySeqOverlaps(readId);
	if (!readId.strand()) readId = readId.rc();
	return _overlapIndex.find(readId).suggestChimeric;
}

std::vector<OverlapRange> 
	OverlapContainer::quickSeqOverlaps(FastaRecord::Id readId, int maxOverlaps, 
									   bool forceLocal)
{
	//bool suggestChimeric;
	const FastaRecord& record = _queryContainer.getRecord(readId);
	return _ovlpDetect.getSeqOverlaps(record, forceLocal, 
									  _divergenceStats, maxOverlaps);
}

std::vector<OverlapRange> 
	OverlapContainer::quickSeqOverlaps(const FastaRecord& record, 
									   int maxOverlaps, bool forceLocal)
{
	return _ovlpDetect.getSeqOverlaps(record, forceLocal, 
									  _divergenceStats, maxOverlaps);
}

const std::vector<OverlapRange>&
	OverlapContainer::lazySeqOverlaps(FastaRecord::Id readId)
{
	bool flipped = !readId.strand();
	if (flipped) readId = readId.rc();
	IndexVecWrapper wrapper;

	//upsert creates default value if it does not exist
	_overlapIndex.upsert(readId, 	
		[&wrapper](IndexVecWrapper& val)
			{wrapper = val;});
	if (wrapper.cached)
	{
		return !flipped ? *wrapper.fwdOverlaps : *wrapper.revOverlaps;
	}

	//otherwise, need to compute overlaps.
	//do it for forward strand to be distinct
	//bool suggestChimeric;
	const bool DEFAULT_LOCAL = false;
	const FastaRecord& record = _queryContainer.getRecord(readId);
	auto overlaps = _ovlpDetect.getSeqOverlaps(record, DEFAULT_LOCAL, 
											   _divergenceStats,
											   _ovlpDetect._maxCurOverlaps);
	overlaps.shrink_to_fit();

	std::vector<OverlapRange> revOverlaps;
	revOverlaps.reserve(overlaps.size());
	for (const auto& ovlp : overlaps) revOverlaps.push_back(ovlp.complement());

	_overlapIndex.update_fn(readId,
		[&wrapper, &overlaps, &revOverlaps, this]
		(IndexVecWrapper& val)
		{
			if (!val.cached)
			{
				_indexSize += overlaps.size();
				*val.fwdOverlaps = std::move(overlaps);
				*val.revOverlaps = std::move(revOverlaps);
				//val.suggestChimeric = suggestChimeric;
				val.cached = true;
			}
			wrapper = val;
		});

	return !flipped ? *wrapper.fwdOverlaps : *wrapper.revOverlaps;
}

void OverlapContainer::ensureTransitivity(bool onlyMaxExt)
{
	Logger::get().debug() << "Computing transitive closure for overlaps";
	
	std::vector<FastaRecord::Id> allSeqs;
	for (const auto& seqIt : _overlapIndex.lock_table()) 
	{
		allSeqs.push_back(seqIt.first);
		allSeqs.push_back(seqIt.first.rc());
	}

	int totalOverlaps = 0;
	for (const auto& seq : allSeqs)
	{
		const auto& curOvlps = this->unsafeSeqOverlaps(seq);
		totalOverlaps += curOvlps.size();
		std::vector<OverlapRange> ovlpsToAdd;
		for (auto& curOvlp : curOvlps)
		{
			auto& extOvlps = this->unsafeSeqOverlaps(curOvlp.extId);

			if (onlyMaxExt)
			{
				bool found = false;
				for (auto& extOvlp : extOvlps)
				{
					if (extOvlp.extId == curOvlp.curId)
					{
						if (curOvlp.score > extOvlp.score)
						{
							extOvlp = curOvlp.reverse();
						}
						found = true;
						break;
					}
				}
				if (!found)
				{
					ovlpsToAdd.push_back(curOvlp.reverse());
				}
			}
			else
			{
				ovlpsToAdd.push_back(curOvlp.reverse());
			}
		}
		for (const auto& ovlp : ovlpsToAdd)
		{
			this->unsafeSeqOverlaps(ovlp.curId).push_back(ovlp);
		}
	}
}


void OverlapContainer::findAllOverlaps()
{
	//Logger::get().info() << "Finding overlaps:";
	std::vector<FastaRecord::Id> allQueries;
	for (const auto& seq : _queryContainer.iterSeqs())
	{
		if (seq.id.strand())
		{
			allQueries.push_back(seq.id);
		}
	}

	std::mutex indexMutex;
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this] (const FastaRecord::Id& seqId)
	{
		this->lazySeqOverlaps(seqId);	//automatically stores overlaps
	};
	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);
	this->ensureTransitivity(false);

	int numOverlaps = 0;
	for (const auto& seqOvlps : _overlapIndex.lock_table()) 
	{
		numOverlaps += seqOvlps.second.fwdOverlaps->size() * 2;
	}
	Logger::get().debug() << "Found " << numOverlaps << " overlaps";

	this->filterOverlaps();

	numOverlaps = 0;
	for (const auto& seqOvlps : _overlapIndex.lock_table()) 
	{
		numOverlaps += seqOvlps.second.fwdOverlaps->size() * 2;
	}
	Logger::get().debug() << "Left " << numOverlaps 
		<< " overlaps after filtering";
}

std::vector<OverlapRange>&
	OverlapContainer::unsafeSeqOverlaps(FastaRecord::Id seqId)
{
		FastaRecord::Id normId = seqId.strand() ? seqId : seqId.rc();
		_overlapIndex.insert(normId);	//ensure it's in the table
		IndexVecWrapper wrapper = _overlapIndex.find(normId);
		return seqId.strand() ? *wrapper.fwdOverlaps : 
								*wrapper.revOverlaps;
}

//TODO: potentially might become non-symmetric after filtering
void OverlapContainer::filterOverlaps()
{
	static const int MAX_ENDS_DIFF = Parameters::get().kmerSize;

	std::vector<FastaRecord::Id> seqIds;
	for (const auto& seq : _queryContainer.iterSeqs())
	{
		seqIds.push_back(seq.id);
	}

	std::function<void(const FastaRecord::Id& seqId)> filterParallel =
	[this] (const FastaRecord::Id& seqId)
	{
		auto& overlaps = this->unsafeSeqOverlaps(seqId);
		
		SetVec<OverlapRange*> overlapSets;
		for (auto& ovlp : overlaps) 
		{
			overlapSets.push_back(new SetNode<OverlapRange*>(&ovlp));
		}
		for (size_t i = 0; i < overlapSets.size(); ++i)
		{
			for (size_t j = 0; j < overlapSets.size(); ++j)
			{
				OverlapRange& ovlpOne = *overlapSets[i]->data;
				OverlapRange& ovlpTwo = *overlapSets[j]->data;

				if (ovlpOne.extId != ovlpTwo.extId) continue;
				int curDiff = ovlpOne.curRange() - ovlpOne.curIntersect(ovlpTwo);
				int extDiff = ovlpOne.extRange() - ovlpOne.extIntersect(ovlpTwo);

				if (curDiff < MAX_ENDS_DIFF && extDiff < MAX_ENDS_DIFF) 
				{
					unionSet(overlapSets[i], overlapSets[j]);
				}
			}
		}
		auto clusters = groupBySet(overlapSets);
		std::vector<OverlapRange> newOvlps;
		for (const auto& cluster : clusters)
		{
			OverlapRange* maxOvlp = nullptr;
			for (auto& ovlp : cluster.second)
			{
				if (!maxOvlp || ovlp->score > maxOvlp->score)
				{
					maxOvlp = ovlp;
				}
			}
			newOvlps.push_back(*maxOvlp);
		}
		newOvlps.shrink_to_fit();
		overlaps = std::move(newOvlps);

		std::sort(overlaps.begin(), overlaps.end(), 
				  [](const OverlapRange& o1, const OverlapRange& o2)
				  {return o1.curBegin < o2.curBegin;});

	};
	processInParallel(seqIds, filterParallel, 
					  Parameters::get().numThreads, false);
}


void OverlapContainer::estimateOverlaperParameters()
{
	Logger::get().debug() << "Estimating k-mer identity bias";

	//const int NEDEED_OVERLAPS = 1000;
	const int MAX_SEQS = 1000;

	std::vector<FastaRecord::Id> readsToCheck;
	for (size_t i = 0; i < MAX_SEQS; ++i) 
	{
		size_t randId = rand() % _queryContainer.iterSeqs().size();
		readsToCheck.push_back(_queryContainer.iterSeqs()[randId].id);
	}

	std::mutex storageMutex;
	std::vector<float> biases;
	std::vector<float> trueDivergence;
	std::function<void(const FastaRecord::Id& seqId)> computeParallel =
	[this, &storageMutex, &biases, &trueDivergence] (const FastaRecord::Id& seqId)
	{
		auto overlaps = this->quickSeqOverlaps(seqId, /*max ovlps*/ 0);
		OverlapRange* maxOvlp = nullptr;
		for (auto& ovlp : overlaps)
		{
			if (!maxOvlp || ovlp.curRange() > maxOvlp->curRange())
			{
				maxOvlp = &ovlp;
			}
		}
		if (maxOvlp)
		{
			std::lock_guard<std::mutex> lock(storageMutex);
			trueDivergence.push_back(maxOvlp->seqDivergence);
		}

		/*for (const auto& ovlp : overlaps)
		{
			float trueDiv = 
				getAlignmentErrEdlib(ovlp, _queryContainer.getSeq(seqId),
									 _ovlpDetect._seqContainer.getSeq(ovlp.extId),
									 0.5f);

			std::lock_guard<std::mutex> lock(storageMutex);
			biases.push_back(trueDiv - ovlp.seqDivergence);
			trueDivergence.push_back(trueDiv);

			std::lock_guard<std::mutex> lock(storageMutex);
			trueDivergence.push_back(ovlp.seqDivergence);
			if (trueDivergence.size() >= NEDEED_OVERLAPS) return;
		}*/
	};
	processInParallel(readsToCheck, computeParallel, 
					  Parameters::get().numThreads, false);

	if (!trueDivergence.empty())
	{
		//_kmerIdyEstimateBias = median(biases);
		//_kmerIdyEstimateBias = 0;
		_meanTrueOvlpDiv = median(trueDivergence);

		//set the parameters and reset statistics
		//_ovlpDetect._estimatorBias = _kmerIdyEstimateBias;
		_divergenceStats.vecSize = 0;
	}
	else
	{
		Logger::get().warning() << "No overlaps found - unable to estimate parameters";
		_meanTrueOvlpDiv = 0.5f;
		//_kmerIdyEstimateBias = 0.0f;
	}

	Logger::get().debug() << "Initial divergence estimate : " << _meanTrueOvlpDiv;
	//Logger::get().debug() << "K-mer estimate bias (true - est): " << _kmerIdyEstimateBias;
}


void OverlapContainer::setDivergenceThreshold(float threshold, bool isRelative)
{

	_ovlpDetect._maxDivergence = (isRelative ? _meanTrueOvlpDiv : 0.0f) + threshold;
	Logger::get().debug() << "Relative threshold: " << "NY"[(size_t)isRelative];
	Logger::get().debug() << "Max divergence threshold set to " 
		<< _ovlpDetect._maxDivergence;
}

void OverlapContainer::overlapDivergenceStats()
{
	this->overlapDivergenceStats(_divergenceStats, 
								 _ovlpDetect._maxDivergence);
}

void OverlapContainer::overlapDivergenceStats(const OvlpDivStats& stats,
											  float divCutoff)
{
	std::vector<float> ovlpDivergence(stats.divVec.begin(),
									  stats.divVec.begin() + stats.vecSize);
	const int HIST_LENGTH = 100;
	const int HIST_HEIGHT = 20;
	const float HIST_MIN = 0;
	const float HIST_MAX = 0.5;
	const float mult = HIST_LENGTH / (HIST_MAX * 100);
	std::vector<int> histogram(HIST_LENGTH, 0);
	for (float d : ovlpDivergence)
	{
		if (HIST_MIN <= d && d < HIST_MAX) 
		{
			++histogram[int(d * mult * 100)];
		}
	}
	int histMax = 1;
	int threshold = divCutoff * mult * 100;
	for (int freq : histogram) histMax = std::max(histMax, freq);

	std::string histString = "\n";
	for (int height = HIST_HEIGHT - 1; height >= 0; --height)
	{
		histString += "    |";
		for (int i = 0; i < HIST_LENGTH; ++i)
		{
			if ((float)histogram[i] / histMax > (float)height / HIST_HEIGHT)
			{
				histString += '*';
			}
			else if (i != threshold)
			{
				histString += ' ';
			}
			else
			{
				histString += '|';
			}
		}
		histString += '\n';
	}
	histString += "    " + std::string(HIST_LENGTH,  '-') + "\n";
	std::string footer(HIST_LENGTH, ' ');
	for (int i = 0; i < 10; ++i)
	{
		size_t startPos = i * HIST_LENGTH / 10;
		auto s = std::to_string(i * 5) + "%";
		for (size_t j = 0; j < s.size(); ++j) footer[j + startPos] = s[j];
	}
	histString += "    " + footer + "\n";

	Logger::get().info() << "Median overlap divergence: " 
		<< quantile(ovlpDivergence, 50); 
	Logger::get().debug() << "Sequence divergence distribution: \n" << histString
		<< "\n    Q25 = " << std::setprecision(2)
		<< quantile(ovlpDivergence, 25) << ", Q50 = " 
		<< quantile(ovlpDivergence, 50)
		<< ", Q75 = " << quantile(ovlpDivergence, 75) << "\n"
		<< std::setprecision(6);
}


void OverlapContainer::buildIntervalTree()
{
	//Logger::get().debug() << "Building interval tree";
	std::vector<FastaRecord::Id> allSeqs;
	for (const auto& seqIt : _overlapIndex.lock_table()) 
	{
		allSeqs.push_back(seqIt.first);
		allSeqs.push_back(seqIt.first.rc());
	}

	for (const auto& seq : allSeqs)
	{
		std::vector<Interval<const OverlapRange*>> intervals;
		auto& overlaps = this->unsafeSeqOverlaps(seq);
		for (const auto& ovlp : overlaps)
		{
			intervals.emplace_back(ovlp.curBegin, ovlp.curEnd, &ovlp);
		}
		_ovlpTree[seq] = IntervalTree<const OverlapRange*>(intervals);
	}
}

std::vector<Interval<const OverlapRange*>> 
	OverlapContainer::getCoveringOverlaps(FastaRecord::Id seqId, 
								  int32_t start, int32_t end) const
{
	return _ovlpTree.at(seqId).findOverlapping(start, end);
}
