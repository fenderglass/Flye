//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "read_aligner.h"
#include "../common/parallel.h"
#include <cmath>
#include <iomanip>

namespace
{
	struct Chain
	{
		Chain(): score(0) {}
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
	//static const int32_t MAX_SEP = Config::get("max_separation");
	static const int32_t MAX_READ_OVLP = 50;
	//static const int32_t KMER_SIZE = Parameters::get().kmerSize;

	std::list<Chain> activeChains;
	for (auto& edgeAlignment : ovlps)
	{
		std::list<Chain> newChains;
		int32_t maxScore = 0;
		Chain* maxChain = nullptr;
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
				//int32_t jumpDiv = abs(readDiff - 
				//					  (graphLeftDiff + graphRightDiff));
				//int32_t gapCost = (jumpDiv > 100) ? 2 * jumpDiv : 0;
				int32_t gapCost = std::max(-readDiff, 0);
				int32_t score = chain.score + nextOvlp.score - gapCost;
				if (score > maxScore)
				{
					maxScore = score;
					maxChain = &chain;
				}
			}
		}
		
		if (maxChain)
		{
			newChains.push_back(*maxChain);
			maxChain->aln.push_back(&edgeAlignment);
			maxChain->score = maxScore;
		}

		activeChains.splice(activeChains.end(), newChains);
		activeChains.push_back(Chain());
		activeChains.back().aln.push_back(&edgeAlignment);
		activeChains.back().score = edgeAlignment.overlap.score;
	}

	//greedily choose non-intersecting set of alignments
	std::vector<GraphAlignment> acceptedAlignments;
	std::vector<Chain> sortedChains(activeChains.begin(), activeChains.end());
	std::sort(sortedChains.begin(), sortedChains.end(),
			  [](const Chain& c1, const Chain& c2)
			  {return c1.score > c2.score;});
	for (auto& chain : sortedChains)
	{
		int32_t alnLen = chain.aln.back()->overlap.curEnd - 
					 	 chain.aln.front()->overlap.curBegin;
		if (alnLen < Parameters::get().minimumOverlap) continue;

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
			if (overlapRate > (int)Config::get("max_separation")) overlaps = true;
		}
		if (!overlaps) 
		{
			acceptedAlignments.emplace_back();
			for (auto& aln : chain.aln) acceptedAlignments.back().push_back(*aln);
		}
	}

	return acceptedAlignments;
}

void ReadAligner::alignReads()
{
	static const int MIN_EDGE_OVLP = (int)Config::get("max_separation");
	static const int EDGE_FLANK = 100;

	//create database
	ska::flat_hash_map<FastaRecord::Id, 
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
	VertexIndex pathsIndex(_graph.edgeSequences(), 
						   (int)Config::get("read_align_kmer_sample"));
	pathsIndex.countKmers(/*min freq*/ 1, /* genome size*/ 0);
	pathsIndex.setRepeatCutoff(/*min freq*/ 1);
	pathsIndex.buildIndex(/*min freq*/ 1);
	OverlapDetector readsOverlapper(_graph.edgeSequences(), pathsIndex, 
									(int)Config::get("maximum_jump"),
									MIN_EDGE_OVLP - EDGE_FLANK,
									/*no overhang*/ 0, /*no max ovlp count*/ 0,
									/*keep alignment*/ false, /*only max*/ false,
									(float)Config::get("read_align_ovlp_divergence"),
									/*bad end adjust*/ 0.0f, 
									/*nucl alignment*/ false);
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

	std::function<void(const FastaRecord::Id&)> alignRead = 
	[this, &indexMutex, &numAligned, &readsOverlaps,
		&idToSegment, &alignedLength, &alignedInFull] 
	(const FastaRecord::Id& seqId)
	{
		auto overlaps = readsOverlaps.quickSeqOverlaps(seqId);
		std::vector<EdgeAlignment> alignments;
		for (auto& ovlp : overlaps)
		{
			//because edges might be as short as max_separation,
			//we set minimum alignment threshold to a bit shorter value.
			//However, apply the actual threshold for longer edges now.
			if (ovlp.extLen < MIN_EDGE_OVLP + EDGE_FLANK ||
				std::min(ovlp.curRange(), ovlp.extRange()) > MIN_EDGE_OVLP)
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

		std::vector<GraphAlignment> complChains(readChains);
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

		if (readChains.empty()) return;

		/////synchronized part
		indexMutex.lock();
		++numAligned;
		if (readChains.size() == 1) ++alignedInFull;
		for (auto& chain : readChains) 
		{
			_readAlignments.push_back(chain);
			alignedLength += chain.back().overlap.curEnd - 
							 chain.front().overlap.curBegin;
		}
		for (auto& chain : complChains)
		{
			_readAlignments.push_back(chain);
		}
		indexMutex.unlock();
		/////
	};

	processInParallel(allQueries, alignRead, 
					  Parameters::get().numThreads, true);

	/*for (auto& aln : _readAlignments)
	{
		if (aln.size() > 1)
		{
			std::string alnStr;
			int switches = 0;
			for (size_t i = 0; i < aln.size() - 1; ++i)
			{
				if (aln[i].segment.end != aln[i + 1].segment.start) ++switches;
			}

			int totalScore = 0;
			int32_t prevGap = 0;
			int32_t prevReadPos = 0;
			for (auto& edge : aln)
			{
				totalScore += edge.overlap.score;
				int32_t nextGap = edge.overlap.extBegin;
				int32_t readGap = edge.overlap.curBegin - prevReadPos;
				if (prevGap > 0)
				{
					int32_t gapCost = (nextGap + prevGap) / 10;
					totalScore -= gapCost;
				}

				alnStr += std::to_string(edge.edge->edgeId.signedId()) + " ("
					+ std::to_string(edge.overlap.curRange()) + ", " 
					+ std::to_string(edge.overlap.seqDivergence) + ", " 
					+ std::to_string(nextGap + prevGap) + ", "
					+ std::to_string(readGap) + ") -> ";

				prevGap = edge.overlap.extLen - edge.overlap.extEnd;
				prevReadPos = edge.overlap.curEnd;
			}
			alnStr.erase(alnStr.size() - 4);
			alnStr += " readLen:" + std::to_string(aln.front().overlap.curLen);
			alnStr += " alnLen:" + std::to_string(aln.back().overlap.curEnd - 
												   aln.front().overlap.curBegin);
			alnStr += " score:" + std::to_string(totalScore);
			FastaRecord::Id readId = aln.front().overlap.curId;
			Logger::get().debug() << "Aln " << _readSeqs.seqName(readId);
			Logger::get().debug() << "\t" << alnStr;
		}
	}*/

	Logger::get().debug() << "Total reads : " << allQueries.size();
	Logger::get().debug() << "Read with aligned parts : " << numAligned;
	Logger::get().debug() << "Aligned in one piece : " << alignedInFull;
	Logger::get().info() << "Aligned read sequence: " << alignedLength << " / " 
		<< totalLength << " (" << (float)alignedLength / totalLength << ")";
	readsOverlaps.overlapDivergenceStats();
}

//updates alignments with respect to the new graph
void ReadAligner::updateAlignments()
{
	std::vector<GraphAlignment> newAlignments;
	for (auto& aln : _readAlignments)
	{
		GraphAlignment curAlignment;
		for (size_t i = 0; i < aln.size() - 1; ++i)
		{
			if (!_graph.hasEdge(aln[i].edge)) continue;

			curAlignment.push_back(aln[i]);
			if (!_graph.hasEdge(aln[i + 1].edge) ||
				aln[i].edge->nodeRight != aln[i + 1].edge->nodeLeft)
			{
				newAlignments.push_back(curAlignment);
				curAlignment.clear();
			}
		}

		if (_graph.hasEdge(aln.back().edge)) curAlignment.push_back(aln.back());
		if (!curAlignment.empty()) newAlignments.push_back(curAlignment);
	}

	_readAlignments = newAlignments;
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
	while(!fin.eof())
	{
		std::string buffer;
		fin >> buffer;
		if (buffer.empty()) continue;

		if (buffer == "Chain")
		{
			if (!curAlignment.empty())
			{
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
			curAlignment.push_back({ovlp, edge});
		}
		else throw std::runtime_error("Error parsing: " + filename);
	}
	if (!curAlignment.empty())
	{
		_readAlignments.push_back(curAlignment);
		curAlignment.clear();
	}

	this->updateAlignments();
}
