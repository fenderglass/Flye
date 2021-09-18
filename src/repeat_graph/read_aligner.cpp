//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "read_aligner.h"
#include "../sequence/alignment.h"
#include "../common/parallel.h"
#include <cmath>
#include <iomanip>
#include <queue>

namespace
{
	struct Chain
	{
		std::vector<const EdgeAlignment*> aln;
		int32_t score;
	};
}

//Give alignments to separate edges for a single read, merges them
//into non-overlapping chains (could be more than one chain per read
//in case of chimera) with maximum score
std::vector<GraphAlignment>
	ReadAligner::chainReadAlignments(const std::vector<EdgeAlignment>& ovlps) const
{
	static const int32_t MAX_JUMP = Config::get("maximum_jump");
	static const int32_t MAX_READ_OVLP = 50;
	static const int32_t MIN_ALN = Parameters::get().minimumOverlap;
	static const int32_t MAX_SEP = (int)Config::get("max_separation");

	std::deque<Chain> activeChains;
	std::deque<Chain> frozenChains;
	for (auto& edgeAlignment : ovlps)
	{
		int32_t maxScore = 0;
		Chain* maxChain = nullptr;
		int numOutdated = 0;

		bool canExtend = edgeAlignment.overlap.extBegin < MAX_JUMP;
		bool canBeExtended = edgeAlignment.overlap.extLen - 
						   	 edgeAlignment.overlap.extEnd < MAX_JUMP;

		if (canExtend)
		{
			for (auto& chain : activeChains)
			{
				const OverlapRange& nextOvlp = edgeAlignment.overlap;
				const OverlapRange& prevOvlp = chain.aln.back()->overlap;

				int32_t readDiff = nextOvlp.curBegin - prevOvlp.curEnd;
				int32_t graphLeftDiff = nextOvlp.extBegin;
				int32_t graphRightDiff = prevOvlp.extLen - prevOvlp.extEnd;

				if (chain.aln.back()->edge->nodeRight == edgeAlignment.edge->nodeLeft &&
					MAX_JUMP > readDiff && readDiff > -MAX_READ_OVLP &&
					graphLeftDiff + graphRightDiff < MAX_JUMP)
				{
					int32_t jumpDiv = abs(readDiff - (graphLeftDiff + graphRightDiff));
					int32_t gapCost = (jumpDiv > 100) ? jumpDiv / 50 : 0;
					int32_t score = chain.score + nextOvlp.score - gapCost;
					if (score > maxScore)
					{
						maxScore = score;
						maxChain = &chain;
					}
				}

				if (readDiff > MAX_JUMP) ++numOutdated;
			}
		}

		//found chain to continue
		if (maxChain)
		{
			activeChains.push_back(*maxChain);		//add a copy of the extended chain
			activeChains.back().aln.push_back(&edgeAlignment);
			activeChains.back().score = maxScore;
		}
		//can't continue, create a new chain
		else
		{
			if (canBeExtended)
			{
				activeChains.push_back({{&edgeAlignment}, 
										 edgeAlignment.overlap.score});
			}
			else
			{
				frozenChains.push_back({{&edgeAlignment}, 
										 edgeAlignment.overlap.score});
			}
		}

		//cleaning up if too much outdated chains
		if (numOutdated > (int)activeChains.size() / 2)
		{
			auto itInsert = activeChains.begin();
			auto itCur = activeChains.begin();
			while(itCur != activeChains.end())
			{
				bool outdated = edgeAlignment.overlap.curBegin - 
								itCur->aln.back()->overlap.curEnd > MAX_JUMP;
				if (outdated)
				{
					frozenChains.push_back(*itCur);
				}
				else
				{
					if (itInsert != itCur) *itInsert = *itCur;
					++itInsert;
				}
				++itCur;
			}
			activeChains.erase(itInsert, activeChains.end());
		}
	}

	activeChains.insert(activeChains.end(), frozenChains.begin(), 
						frozenChains.end());	
	std::sort(activeChains.begin(), activeChains.end(),
			  [](const Chain& c1, const Chain& c2)
			  {return c1.score > c2.score;});

	//greedily choose non-intersecting set of alignments
	std::vector<GraphAlignment> acceptedAlignments;
	for (auto& chain : activeChains)
	{
		int32_t alnLen = chain.aln.back()->overlap.curEnd - 
					 	 chain.aln.front()->overlap.curBegin;
		if (alnLen < MIN_ALN) continue;

		//check if it overlaps with other accepted chains
		bool overlaps = false;
		for (auto& existAln : acceptedAlignments)
		{
			int32_t existStart = existAln.front().overlap.curBegin;
			int32_t existEnd = existAln.back().overlap.curEnd;
			int32_t curStart = chain.aln.front()->overlap.curBegin;
			int32_t curEnd = chain.aln.back()->overlap.curEnd;

			int32_t overlapRate = std::min(curEnd, existEnd) - 
									std::max(curStart, existStart);
			if (overlapRate > MAX_SEP) overlaps = true;
		}
		if (!overlaps) 
		{
			acceptedAlignments.emplace_back();
			for (auto& aln : chain.aln) acceptedAlignments.back().push_back(*aln);
		}
	}

	return acceptedAlignments;
}

void ReadAligner::alignReads(const std::string& kmerWhitelistPath)
{
	static const int SMALL_ALN = 100;
	static const int BIG_ALN = 500;
	static const int LONG_EDGE = 900;

	static const float MAX_DIVERGENCE = Config::get("read_align_ovlp_divergence");

	//create database
	std::unordered_map<FastaRecord::Id, 
					   std::pair<GraphEdge*, EdgeSequence>> idToSegment;
	for (auto& edge : _graph.iterEdges())
	{
		for (auto& segment : edge->seqSegments)
		{
			idToSegment[segment.edgeSeqId] = {edge, segment};
			idToSegment[segment.edgeSeqId.rc()] = {_graph.complementEdge(edge), 
										   		   segment.complement()};
		}
	}

	//index it and align reads
	VertexIndex pathsIndex(_graph.edgeSequences());
	bool useMinimizers = Config::get("use_minimizers");
	int minWnd = useMinimizers ? Config::get("minimizer_window") : 1;
	pathsIndex.buildIndexMinimizers(/*min freq*/ 1, minWnd);

	if (!kmerWhitelistPath.empty())
	{
		pathsIndex.loadKmerWhitelist(kmerWhitelistPath);
	}

	//pathsIndex.countKmers(/*min freq*/ 1, /* genome size*/ 0);
	//pathsIndex.buildIndex(/*min freq*/ 1);
	OverlapDetector readsOverlapper(_graph.edgeSequences(), pathsIndex, 
									(int)Config::get("maximum_jump"), SMALL_ALN,
									/*no overhang*/ 0, /*keep alignment*/ false, 
									/*only max*/ false, /*no max divergence*/ 1.0f,
									/*nucl alignment*/ false,
									/*partition bad map*/ false,
								    (bool)Config::get("hpc_scoring_on"));
	OverlapContainer readsOverlaps(readsOverlapper, _readSeqs);

	std::vector<FastaRecord::Id> allQueries;
	int64_t totalLength = 0;
	for (auto& read : _readSeqs.iterSeqs())
	{
		if (!read.id.strand()) continue;
		if (read.sequence.length() > (size_t)Parameters::get().minimumOverlap)
		{
			totalLength += read.sequence.length();
			allQueries.push_back(read.id);
		}
	}
	std::mutex indexMutex;
	int numAligned = 0;
	int alignedInFull = 0;
	int64_t alignedLength = 0;
	OvlpDivStats divergenceStats;

	std::function<void(const FastaRecord::Id&)> alignRead = 
	[this, &indexMutex, &numAligned, &readsOverlaps,
		&idToSegment, &alignedLength, &alignedInFull, &divergenceStats] 
	(const FastaRecord::Id& seqId)
	{
		auto overlaps = readsOverlaps.quickSeqOverlaps(seqId);
		std::vector<EdgeAlignment> alignments;
		for (auto& ovlp : overlaps)
		{
			//because edges might be as short as max_separation,
			//we set minimum alignment threshold to a bit shorter value.
			//However, apply the actual threshold for longer edges now.
			if (ovlp.extLen < LONG_EDGE ||
				std::min(ovlp.curRange(), ovlp.extRange()) > BIG_ALN)
			{
				//alignments.push_back({ovlp, idToSegment[ovlp.extId].first,
				//					  idToSegment[ovlp.extId].second});
				alignments.push_back({ovlp, idToSegment[ovlp.extId].first});
			}

		}
		std::sort(alignments.begin(), alignments.end(),
		  [](const EdgeAlignment& e1, const EdgeAlignment& e2)
			{return e1.overlap.curBegin < e2.overlap.curBegin;});
		auto readChains = this->chainReadAlignments(alignments);

		//check divergence once the chain is formed
		std::vector<GraphAlignment> goodChains;
		for (auto& chain : readChains)
		{
			float chainDivergence = 
				this->getChainBaseDivergence(chain, (bool)Config::get("reads_base_alignment"));
			divergenceStats.add(chainDivergence);
			if (chainDivergence < MAX_DIVERGENCE)
			{
				goodChains.push_back(chain);
			}
		}

		std::vector<GraphAlignment> complChains(goodChains);
		for (auto& chain : complChains)
		{
			for (auto& aln : chain)
			{
				aln.edge = _graph.complementEdge(aln.edge);
				//aln.segment = aln.segment.complement();
				aln.overlap = aln.overlap.complement();
			}
			std::reverse(chain.begin(), chain.end());
		}

		if (goodChains.empty()) return;

		/////synchronized part
		indexMutex.lock();
		++numAligned;
		if (goodChains.size() == 1) ++alignedInFull;
		for (auto& chain : goodChains) 
		{
			chain.shrink_to_fit();
			_readAlignments.push_back(chain);
			alignedLength += chain.back().overlap.curEnd - 
							 chain.front().overlap.curBegin;
		}
		for (auto& chain : complChains)
		{
			chain.shrink_to_fit();
			_readAlignments.push_back(chain);
		}
		indexMutex.unlock();
		/////
	};

	processInParallel(allQueries, alignRead, 
					  Parameters::get().numThreads, true);

	Logger::get().debug() << "Total reads : " << allQueries.size();
	Logger::get().debug() << "Read with aligned parts : " << numAligned;
	Logger::get().debug() << "Aligned in one piece : " << alignedInFull;
	Logger::get().info() << "Aligned read sequence: " << alignedLength << " / " 
		<< totalLength << " (" << (float)alignedLength / totalLength << ")";
	readsOverlaps.overlapDivergenceStats(divergenceStats, MAX_DIVERGENCE);
}

//updates alignments with respect to the new graph
void ReadAligner::updateAlignments()
{
	auto isValidAlignment = [this](const GraphAlignment& aln)
	{
		for (size_t i = 0; i < aln.size() - 1; ++i)
		{
			if (!_graph.getEdge(aln[i].edge->edgeId) ||
				!_graph.getEdge(aln[i + 1].edge->edgeId)) return false;

			if (aln[i].edge->nodeRight != aln[i + 1].edge->nodeLeft) return false;
		}
		return true;
	};

	std::vector<GraphAlignment> newlyAdded;
	auto splitAlignment = [&newlyAdded, this](const GraphAlignment& aln)
	{
		GraphAlignment curAlignment;
		for (size_t i = 0; i < aln.size() - 1; ++i)
		{
			if (!_graph.getEdge(aln[i].edge->edgeId)) continue;

			curAlignment.push_back(aln[i]);
			if (!_graph.getEdge(aln[i + 1].edge->edgeId) ||
				aln[i].edge->nodeRight != aln[i + 1].edge->nodeLeft)
			{
				newlyAdded.push_back(curAlignment);
				curAlignment.clear();
			}
		}

		if (_graph.getEdge(aln.back().edge->edgeId)) curAlignment.push_back(aln.back());
		if (!curAlignment.empty()) newlyAdded.push_back(curAlignment);
	};

	size_t insertIdx = 0;
	for (size_t i = 0; i < _readAlignments.size(); ++i)
	{
		if (isValidAlignment(_readAlignments[i]))
		{
			if (i != insertIdx)
			{
				_readAlignments[insertIdx] = std::move(_readAlignments[i]);
			}
			++insertIdx;
		}
		else
		{
			splitAlignment(_readAlignments[i]);
		}
	}
	_readAlignments.erase(_readAlignments.begin() + insertIdx, _readAlignments.end());
	_readAlignments.reserve(_readAlignments.size() + newlyAdded.size());
	for (auto& aln : newlyAdded)
	{
		_readAlignments.push_back(std::move(aln));
	}
}

void ReadAligner::storeAlignments(const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout)
	{
		throw std::runtime_error("Can't open "  + filename);
	}

	for (auto& chain : _readAlignments)
	{
		fout << "Chain\n";
		for (auto& aln : chain)
		{
			fout << "\tAln\t" << aln.edge->edgeId << "\t";
			aln.overlap.dump(fout, _readSeqs, _graph.edgeSequences());
			fout << "\n";
		}
	}
}

void ReadAligner::loadAlignments(const std::string& filename)
{
	std::ifstream fin(filename);
	if (!fin)
	{
		throw std::runtime_error("Can't open "  + filename);
	}

	GraphAlignment curAlignment;
	while(true)
	{
		std::string buffer;
		fin >> buffer;
		if (fin.eof()) break;
		if (!fin.good()) throw std::runtime_error("Error parsing: " + filename);

		if (buffer == "Chain")
		{
			if (!curAlignment.empty())
			{
				curAlignment.shrink_to_fit();
				_readAlignments.push_back(curAlignment);
				curAlignment.clear();
			}
		}
		else if (buffer == "Aln")
		{
			OverlapRange ovlp;
			size_t edgeId = 0;
			fin >> edgeId;
			ovlp.load(fin, _readSeqs, _graph.edgeSequences());
			GraphEdge* edge = _graph.getEdge(FastaRecord::Id(edgeId));
			if (edge) 
			{
				//sometimes alignment might contain edges that were 
				//removed from the graph (for example, after Trestle).
				//so, we check if the edge exists
				curAlignment.push_back({ovlp, edge});
			}
		}
		else throw std::runtime_error("Error parsing: " + filename);
	}
	if (!curAlignment.empty())
	{
		curAlignment.shrink_to_fit();
		_readAlignments.push_back(curAlignment);
		curAlignment.clear();
	}

	this->updateAlignments();
}

ReadAligner::AlnIndex ReadAligner::makeAlignmentIndex()
{
	AlnIndex alnIndex;
	for (auto& aln : this->getAlignments())
	{
		if (aln.size() > 1)
		{
			std::unordered_set<GraphEdge*> uniqueEdges;
			for (auto& edgeAln : aln)
			{
				uniqueEdges.insert(edgeAln.edge);
			}
			for (GraphEdge* edge : uniqueEdges) alnIndex[edge].push_back(aln);
		}
	}
	return alnIndex;
}

float ReadAligner::getChainBaseDivergence(const GraphAlignment& chain, bool realign)
{
	static const float MAX_DIVERGENCE = Config::get("read_align_ovlp_divergence");
	static const bool USE_HPC = (bool)Config::get("hpc_scoring_on");

	float sumMatched = 0;
	int alnLen = 0;
	for (auto& aln : chain)
	{
		float ovlpDivergence = aln.overlap.seqDivergence;
		if (realign)
		{
			ovlpDivergence = 
				getAlignmentErrEdlib(aln.overlap, _readSeqs.getSeq(aln.overlap.curId), 
									 _graph.edgeSequences().getSeq(aln.overlap.extId),
									 MAX_DIVERGENCE, USE_HPC);
		}

		sumMatched += aln.overlap.curRange() * (1 - ovlpDivergence);
		alnLen += aln.overlap.curRange();
	}
	float chainDivergence = 1 - (float)sumMatched / alnLen;
	
	return chainDivergence;
}


ReadAligner::ConnIndex ReadAligner::getEdgeConnectivity() const
{
	ConnIndex connections;
	for (auto& aln : _readAlignments)
	{
		for (size_t i = 0; i < aln.size() - 1; ++i)
		{
			++connections[aln[i].edge][aln[i + 1].edge];
		}
	}
	return connections;
}
