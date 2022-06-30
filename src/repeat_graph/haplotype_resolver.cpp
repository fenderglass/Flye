#include "haplotype_resolver.h"
#include "graph_processing.h"
#include <queue>
#include <set>

//This function collapses simple bubbles caused by
//alternative haplotypes / strains. They are defined as follows:
//1. Structure: 1 input, 2 branches, 1 output: -<>-
//2. Size of each branch is shorter than MAX_BUBBLE_LEN below
//3. Each branch is shorter than both entrace and exits. We need this to
//   distinguish from the case of two repeats of multiplicity 2
//Note that we are not using any global coverage assumptions here.
int HaplotypeResolver::findHeterozygousBulges()
{
	//const float MAX_COV_VAR = 1.5;
	const int MAX_BUBBLE_LEN = Config::get("max_bubble_length");

	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	std::unordered_map<GraphEdge*, UnbranchingPath*> pathIndex;
	for (UnbranchingPath& path : unbranchingPaths)
	{
		for (GraphEdge* edge : path.path) pathIndex[edge] = &path;
	}

	std::unordered_set<FastaRecord::Id> usedPaths;
	int numMasked = 0;
	for (auto& path : unbranchingPaths)
	{
		if (path.isLooped()) continue;
		
		if (path.nodeLeft()->inEdges.size() != 1 ||
			path.nodeLeft()->outEdges.size() != 2 ||
			path.nodeRight()->outEdges.size() != 1 ||
			path.nodeRight()->inEdges.size() != 2) continue;

		std::vector<UnbranchingPath*> twoPaths;
		for (GraphEdge* edge : path.nodeLeft()->outEdges)
		{
			if (pathIndex[edge]->nodeRight() == path.nodeRight())
			{
				twoPaths.push_back(pathIndex[edge]);
			}
		}
		if (twoPaths.size() != 2) continue;
		if (twoPaths[0]->id == twoPaths[1]->id.rc()) continue;
		if (usedPaths.count(twoPaths[0]->id) || 
			usedPaths.count(twoPaths[1]->id)) continue;

		UnbranchingPath* entrancePath = pathIndex[path.nodeLeft()->inEdges[0]];
		UnbranchingPath* exitPath = pathIndex[path.nodeRight()->outEdges[0]];

		//if (entrancePath->id == exitPath->id ||
		if (entrancePath->id == exitPath->id.rc()) continue;

		//sanity check for maximum bubble size
		if (std::max(twoPaths[0]->length, twoPaths[1]->length) > 
			MAX_BUBBLE_LEN) continue;

		//coverage requirement: sum over two branches roughly equals to
		//exit and entrance coverage or less
		//float covSum = twoPaths[0]->meanCoverage + twoPaths[1]->meanCoverage;
		//if (covSum > std::min(entrancePath->meanCoverage * MAX_COV_VAR,
		//					  exitPath->meanCoverage * MAX_COV_VAR)) continue;

		//require bubble branches to be shorter than entrance or exit,
		//to distinguish from the case of two consecutive repeats
		//of multiplicity 2
		//if (std::max(twoPaths[0]->length, twoPaths[1]->length) >
		//	std::max(entrancePath->length, exitPath->length)) continue;
		
		for (size_t i = 0; i < 2; ++i)
		{
			usedPaths.insert(twoPaths[i]->id);
			usedPaths.insert(twoPaths[i]->id.rc());
		}

		if (twoPaths[0]->meanCoverage > twoPaths[1]->meanCoverage)
		{
			std::swap(twoPaths[0], twoPaths[1]);
		}

		if (!twoPaths[0]->path.front()->altHaplotype ||
			!twoPaths[1]->path.front()->altHaplotype) ++numMasked;

		for (size_t i = 0; i < 2; ++i)
		{
			for (auto& edge : twoPaths[i]->path)
			{
				edge->altHaplotype = true;
				edge->altGroupId = _nextAltGroupId;
				_graph.complementEdge(edge)->altHaplotype = true;
				_graph.complementEdge(edge)->altGroupId = _nextAltGroupId + 1;
			}
		}
		_nextAltGroupId += 2;

		GraphEdge* inEdge = entrancePath->path.back();
		GraphEdge* outEdge = exitPath->path.front();

		Logger::get().debug() << "Regular bubble: " << inEdge->edgeDescrLong()
			<< " " << outEdge->edgeDescrLong();
		Logger::get().debug() << "\tInternal: " << twoPaths[0]->edgesStrLong() << " " << twoPaths[1]->edgesStrLong();

		if (inEdge->rightLink || outEdge->leftLink) continue;


		//link edges
		_graph.linkEdges(inEdge, outEdge);
		_graph.linkEdges(_graph.complementEdge(outEdge),
						 _graph.complementEdge(inEdge));

		//bridging sequence
		DnaSequence pathSeq = this->pathSequence(twoPaths[0]->path);
		_bridgingSeqs[std::make_pair(inEdge, outEdge)] = pathSeq;
		_bridgingSeqs[std::make_pair(_graph.complementEdge(outEdge), 
									 _graph.complementEdge(inEdge))] = 
												pathSeq.complement();

	}

	Logger::get().debug() << "[SIMPL] Masked " << numMasked
		<< " simple bubbles";
	return numMasked;
}

//This function collapses simple loops:
//1. One loop edge with one entrance and one exit
//2. Loop length is shorter than lengths of entrance/exit
//3. Loop coverage is roughly equal or less than coverage of entrance/exit
int HaplotypeResolver::findHeterozygousLoops()
{
	const float COV_MULT = (float)Config::get("loop_coverage_rate");
	const int MAX_LOOP_LEN = Config::get("max_bubble_length");

	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	std::unordered_set<FastaRecord::Id> toUnroll;
	std::unordered_set<FastaRecord::Id> toRemove;
	int numMasked = 0;
	for (auto& loop : unbranchingPaths)
	{
		if (!loop.id.strand()) continue;
		if (!loop.isLooped()) continue;
		if (loop.path.front()->selfComplement) continue;
		if (loop.length > MAX_LOOP_LEN) continue;

		GraphNode* node = loop.nodeLeft();
		if (node->inEdges.size() != 2 ||
			node->outEdges.size() != 2) continue;

		UnbranchingPath* entrancePath = nullptr;
		UnbranchingPath* exitPath = nullptr;
		for (auto& cand : unbranchingPaths)
		{
			if (cand.nodeRight() == node &&
				loop.id != cand.id) entrancePath = &cand;
			if (cand.nodeLeft() == node &&
				loop.id != cand.id) exitPath = &cand;
		}

		if (entrancePath->isLooped()) continue;
		if (entrancePath->id == exitPath->id.rc()) continue;

		//loop coverage should be roughly equal or less
		if (loop.meanCoverage > 
				COV_MULT * std::min(entrancePath->meanCoverage, 
									entrancePath->meanCoverage)) continue;

		//loop should not be longer than other branches
		if (loop.length > std::max(entrancePath->length, 
								   exitPath->length)) continue;

		if (!loop.path.front()->altHaplotype) ++numMasked;
		for (auto& edge : loop.path)
		{
			edge->altHaplotype = true;
			edge->altGroupId = _nextAltGroupId;
			_graph.complementEdge(edge)->altHaplotype = true;
			_graph.complementEdge(edge)->altGroupId = _nextAltGroupId + 1;
		}
		_nextAltGroupId += 2;

		GraphEdge* inEdge = entrancePath->path.back();
		GraphEdge* outEdge = exitPath->path.front();
		if (inEdge->rightLink || outEdge->leftLink) continue;

		Logger::get().debug() << "Bubble-loop: " << inEdge->edgeDescrLong()
			<< " " << outEdge->edgeDescrLong();
		Logger::get().debug() << "\tInternal: " <<  loop.edgesStrLong();

		//links
		_graph.linkEdges(inEdge, outEdge);
		_graph.linkEdges(_graph.complementEdge(outEdge),
						 _graph.complementEdge(inEdge));

		//either remove or unroll loop, depending on the coverage
		int32_t mainPathCoverage = (entrancePath->meanCoverage + exitPath->meanCoverage) / 2;
		int32_t LOOP_COV_RATE = 4;	//2 times lower than half coverage (in case of haploid indel)
		if (loop.meanCoverage < mainPathCoverage / LOOP_COV_RATE)
		{
			_bridgingSeqs[std::make_pair(inEdge, outEdge)] = DnaSequence("A");
			_bridgingSeqs[std::make_pair(_graph.complementEdge(outEdge), 
							  			  _graph.complementEdge(inEdge))] =
															DnaSequence("A");

		}
		else
		{
			DnaSequence seq = this->pathSequence(loop.path);
			_bridgingSeqs[std::make_pair(inEdge, outEdge)] = seq;
			_bridgingSeqs[std::make_pair(_graph.complementEdge(outEdge), 
							  			  _graph.complementEdge(inEdge))] =
														seq.complement();
		}
	}

	Logger::get().debug() << "[SIMPL] Masked " << numMasked << " heterozygous loops"; 
	return numMasked;
}

HaplotypeResolver::VariantPaths 
	HaplotypeResolver::findVariantSegment(GraphEdge* startEdge,
										  const std::vector<GraphAlignment>& alingnments,
										  const std::unordered_set<GraphEdge*>& loopedEdges)
{
	//first, extract alnignment paths starting from
	//the current edge and sort them from longest to shortest
	std::vector<GraphAlignment> outPaths;
	for (auto& aln : alingnments)
	{
		for (size_t i = 0; i < aln.size(); ++i)
		{
			//if (aln[i].edge == startEdge)
			if (aln[i].edge == startEdge && i + 1 < aln.size())
			{
				outPaths.emplace_back(GraphAlignment(aln.begin() + i, 
													 aln.end()));
				break;
			}
		}
	}
	if (outPaths.empty()) return VariantPaths();

	std::sort(outPaths.begin(), outPaths.end(),
			  [](const GraphAlignment& a1, const GraphAlignment& a2)
			  {return a1.back().overlap.curEnd - a1.front().overlap.curEnd >
					  a2.back().overlap.curEnd - a2.front().overlap.curEnd;});

	//Logger::get().debug() << "Haplo paths " 
	//	<< startEdge->edgeId.signedId() << " " << outPaths.size();
	/*for (auto& aln : outPaths)
	{
		std::string pathStr;
		for (size_t i = 0; i < aln.size(); ++i)
		{
			pathStr += std::to_string(aln[i].edge->edgeId.signedId()) + " -> ";
		}
		Logger::get().debug() << "\tPath: " << pathStr;
	}*/

	//now group the path by containmnent. For each group we'll have 
	//a longest "reference" path.

	const int MIN_SCORE = 2;
	std::vector<PathWithScore> pathGroups;
	for (auto& trgPath: outPaths)
	{
		bool newPath = true;
		for (auto& referencePath : pathGroups)
		{
			bool contained = true;
			for (size_t i = 0; i < std::min(trgPath.size(), 
											referencePath.path.size()); ++i)
			{
				if (trgPath[i].edge != referencePath.path[i].edge)
				{
					contained = false;
					break;
				}
			}
			if (contained)
			{
				newPath = false;
				++referencePath.score;
				break;
			}
		}
		if (newPath)
		{
			pathGroups.push_back({trgPath, 1});
		}
	}

	pathGroups.erase(std::remove_if(pathGroups.begin(), pathGroups.end(),
					 [](PathWithScore& p)
					 {return p.score < MIN_SCORE;}), pathGroups.end());

	/*for (auto& aln : pathGroups)
	{
		std::string pathStr;
		for (size_t i = 0; i < aln.path.size(); ++i)
		{
			pathStr += std::to_string(aln.path[i].edge->edgeId.signedId()) + " -> ";
		}
		Logger::get().debug() << "\tGroup: " << pathStr << aln.score;
	}*/

	if (pathGroups.size() < 2) return VariantPaths();

	//mark edges that appear more than once as repeats
	std::unordered_set<GraphEdge*> repeats;
	for (size_t groupId = 0; groupId < pathGroups.size(); ++groupId)
	{
		std::unordered_set<GraphEdge*> seen;
		for (size_t i = 0; i < pathGroups[groupId].path.size(); ++i)
		{
			if (seen.count(pathGroups[groupId].path[i].edge))
			{
				repeats.insert(pathGroups[groupId].path[i].edge);
			}
			seen.insert(pathGroups[groupId].path[i].edge);
		}
	}

	//now, set the longest path as reference, and find
	//edges where other groups coverge with the reference
	PathWithScore& refPath = pathGroups.front();
	std::unordered_set<GraphEdge*> convergenceEdges;
	for (size_t i = 0; i < refPath.path.size(); ++i)
	{
		if (!loopedEdges.count(refPath.path[i].edge) &&
			!repeats.count(refPath.path[i].edge))
		{
			convergenceEdges.insert(refPath.path[i].edge);
		}
	}
	for (size_t groupId = 1; groupId < pathGroups.size(); ++groupId)
	{
		std::unordered_set<GraphEdge*> newSet;
		for (size_t i = 0; i < pathGroups[groupId].path.size(); ++i)
		{
			if (convergenceEdges.count(pathGroups[groupId].path[i].edge))
			{
				newSet.insert(pathGroups[groupId].path[i].edge);
			}
		}
		convergenceEdges = newSet;
	}

	//get the bubble start (paths might be convergent for a bit)
	size_t bubbleStartId = 0;
	for (;;)
	{
		bool agreement = true;
		for (size_t groupId = 1; groupId < pathGroups.size(); ++groupId)
		{
			if (bubbleStartId + 1 >= pathGroups[groupId].path.size() ||
				!convergenceEdges.count(pathGroups[0].path[bubbleStartId + 1].edge) ||
					(pathGroups[groupId].path[bubbleStartId + 1].edge !=
					 pathGroups[0].path[bubbleStartId + 1].edge))
			{
				agreement = false;
				break;
			}
		}
		if (!agreement) break;
		++bubbleStartId;
	}
	if (!convergenceEdges.count(refPath.path[bubbleStartId].edge)) return VariantPaths();

	//get the bubble end
	bool foundEnd = false;
	size_t bubbleEndId = bubbleStartId + 1;
	for (; bubbleEndId < refPath.path.size(); ++bubbleEndId)
	{
		if (convergenceEdges.count(refPath.path[bubbleEndId].edge))
		{
			foundEnd = true;
			break;
		}
	}
	if (!foundEnd) return VariantPaths();

	//shorten all branches accordingly
	std::vector<PathWithScore> bubbleBranches;
	for (size_t groupId = 0; groupId < pathGroups.size(); ++groupId)
	{
		size_t groupStart = 0;
		size_t groupEnd = 0;
		for (size_t i = 0; i < pathGroups[groupId].path.size(); ++i)
		{
			if (pathGroups[groupId].path[i].edge == 
				refPath.path[bubbleStartId].edge) groupStart = i;

			if (pathGroups[groupId].path[i].edge == 
				refPath.path[bubbleEndId].edge) groupEnd = i;
		}
		GraphAlignment newPath(pathGroups[groupId].path.begin() + groupStart,
							   pathGroups[groupId].path.begin() + groupEnd + 1);
		PathWithScore newBranch = {newPath, pathGroups[groupId].score};

		bool duplicate = false;
		for (size_t branchId = 0; branchId < bubbleBranches.size(); ++branchId)
		{
			if (newBranch.path.size() != 
				bubbleBranches[branchId].path.size()) continue;
			if (std::equal(newBranch.path.begin(), newBranch.path.end(),
						   bubbleBranches[branchId].path.begin(),
						   [](EdgeAlignment& a1, EdgeAlignment& a2)
						   {return a1.edge == a2.edge;}))
			{
				duplicate = true;
				bubbleBranches[branchId].score += newBranch.score;
			}
		}
		if (!duplicate) bubbleBranches.push_back(newBranch);
	}
	if (bubbleBranches.size() < 2) return VariantPaths();

	/*for (auto& aln : bubbleBranches)
	{
		std::string pathStr;
		for (size_t i = 0; i < aln.path.size(); ++i)
		{
			pathStr += std::to_string(aln.path[i].edge->edgeId.signedId()) + " -> ";
		}
		Logger::get().debug() << "\tBranch: " << pathStr << aln.score;
	}*/

	VariantPaths vp;
	vp.startEdge = refPath.path[bubbleStartId].edge;
	vp.endEdge = refPath.path[bubbleEndId].edge;
	vp.altPaths = bubbleBranches;

	//get the bridgin read sequence
	std::vector<GraphAlignment> bridgingReads;
	for (auto& aln : alingnments)
	{
		int startPos = -1;
		int endPos = -1;
		for (size_t i = 0; i < aln.size(); ++i)
		{
			if (aln[i].edge == startEdge) startPos = i;
			if (startPos != -1 && aln[i].edge == vp.endEdge)
			{
				endPos = i;
				break;
			}
		}
		if (startPos != -1 && endPos != -1)
		{
			bridgingReads.emplace_back(aln.begin() + startPos, 
									   aln.begin() + endPos + 1);
		}
	}
	if (bridgingReads.empty()) throw std::runtime_error("No bridging reads!");
	std::sort(bridgingReads.begin(), bridgingReads.end(),
			  [](const GraphAlignment& a1, const GraphAlignment& a2)
			  {return a1.back().overlap.curBegin - a1.front().overlap.curEnd <
			  		  a2.back().overlap.curBegin - a2.front().overlap.curEnd;});

	size_t medianRead = bridgingReads.size() / 2;
	auto readId = bridgingReads[medianRead][0].overlap.curId;
	int32_t readStart = bridgingReads[medianRead][0].overlap.curEnd;
	int32_t readEnd = bridgingReads[medianRead].back().overlap.curBegin;
	const int MAGIC_100 = 100;
	readEnd = std::max(readStart + MAGIC_100 - 1, readEnd);	
	vp.bridgingSequence = _readSeqs.getSeq(readId)
									.substr(readStart, readEnd - readStart);

	return vp;
}

//this function reveals complex heterogenities on the graph
//(more than just two alternative branches) using read-paths
int HaplotypeResolver::findRoundabouts()
{
	auto alnIndex = _aligner.makeAlignmentIndex();

	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	std::unordered_set<GraphEdge*> loopedEdges;
	for (auto& path : unbranchingPaths)
	{
		if (path.isLooped())
		{
			loopedEdges.insert(path.path.begin(), path.path.end());
		}
	}

	std::unordered_set<GraphEdge*> usedEdges;
	std::vector<VariantPaths> foundVariants;
	for (auto& startPath: unbranchingPaths)
	{
		GraphEdge* startEdge = startPath.path.back();
		//if (startEdge->nodeRight->outEdges.size() < 2) continue;
		if (loopedEdges.count(startEdge)) continue;
		if (usedEdges.count(startEdge)) continue;
		
		auto varSeg = this->findVariantSegment(startEdge, alnIndex[startEdge], 
											   loopedEdges);
		if (varSeg.startEdge && varSeg.endEdge &&
			varSeg.startEdge != _graph.complementEdge(varSeg.endEdge))
		{
			auto revSeg = 
				this->findVariantSegment(_graph.complementEdge(varSeg.endEdge), 
										 alnIndex[_graph.complementEdge(varSeg.endEdge)], 
										 loopedEdges);
			if (revSeg.endEdge == _graph.complementEdge(varSeg.startEdge))
			{
				foundVariants.push_back(varSeg);
				usedEdges.insert(revSeg.startEdge);
			}
		}
	}

	int foundNew = 0;
	for (auto& varSegment : foundVariants)
	{
		bool newVariant = true;
		for (auto& branch : varSegment.altPaths)
		{
			for (size_t i = 1; i < branch.path.size() - 1; ++i)
			{
				if (branch.path[i].edge->altHaplotype) newVariant = false;
			}
		}
		if (newVariant)
		{
			++foundNew;
			Logger::get().debug() << "Roundabout: " 
				<< varSegment.startEdge->edgeDescrLong()
				<< " : " << varSegment.endEdge->edgeDescrLong();

			std::string pathStr;
			for (auto& branch : varSegment.altPaths)
			{
				for (size_t i = 1; i < branch.path.size() - 1; ++i)
				{
					pathStr += branch.path[i].edge->edgeDescrLong() + " ";
				}
			}
			Logger::get().debug() << "\tInternal: " << pathStr;
		}

		for (auto& branch : varSegment.altPaths)
		{
			for (size_t i = 1; i < branch.path.size() - 1; ++i)
			{
				branch.path[i].edge->altHaplotype = true;
				branch.path[i].edge->altGroupId = _nextAltGroupId;
				_graph.complementEdge(branch.path[i].edge)->altHaplotype = true;
				_graph.complementEdge(branch.path[i].edge)->altGroupId = _nextAltGroupId + 1;
			}
		}
		_nextAltGroupId += 2;

		if (varSegment.startEdge->rightLink || 
			varSegment.endEdge->leftLink) continue;

		//add links
		_graph.linkEdges(varSegment.startEdge, varSegment.endEdge);
		_graph.linkEdges(_graph.complementEdge(varSegment.endEdge), 
						 _graph.complementEdge(varSegment.startEdge));

		auto fwdPair = std::make_pair(varSegment.startEdge, varSegment.endEdge);
		auto revPair = std::make_pair(_graph.complementEdge(varSegment.endEdge), 
									  _graph.complementEdge(varSegment.startEdge));
		_bridgingSeqs[fwdPair] = varSegment.bridgingSequence;
		_bridgingSeqs[revPair] = varSegment.bridgingSequence.complement();
	}

	Logger::get().debug() << "[SIMPL] Masked " << foundNew << " roundabouts";
	return foundVariants.size();
}

void HaplotypeResolver::collapseHaplotypes()
{
	int numBridged = 0;
	std::unordered_set<GraphEdge*> separatedEdges;
	for (auto& inEdge : _graph.iterEdges())
	{
		if (!inEdge->rightLink) continue;
		if (separatedEdges.count(inEdge)) continue;

		GraphEdge* outEdge = inEdge->rightLink;
		if (!_graph.getEdge(outEdge->edgeId))
		{
			Logger::get().warning() << "Missing linked edge";
			continue;
		}
		if (outEdge->leftLink != inEdge)
		{
			Logger::get().warning() << "Broken link";
			continue; 
		}

		if (!_bridgingSeqs.count(std::make_pair(inEdge, outEdge)))
		{
			Logger::get().warning() << "No bridging path!";
			continue;
		}

		++numBridged;
		separatedEdges.insert(_graph.complementEdge(outEdge));

		DnaSequence& insertSeq = _bridgingSeqs[std::make_pair(inEdge, outEdge)];
		if (insertSeq.length() < 10)	//arbirtary small number
		{
			this->separeteAdjacentEdges(inEdge, outEdge);
			this->separeteAdjacentEdges(_graph.complementEdge(outEdge),
										_graph.complementEdge(inEdge));
		}
		else
		{

			FastaRecord::Id edgeId = _graph.newEdgeId();
			std::stringstream ss;
				ss << "edge_" << edgeId.signedId() << "_haplotype";
			EdgeSequence edgeSeq = 
				_graph.addEdgeSequence(insertSeq, 0, insertSeq.length(), ss.str());

			this->separateDistantEdges(inEdge, outEdge, edgeSeq, edgeId);
			this->separateDistantEdges(_graph.complementEdge(outEdge),
									   _graph.complementEdge(inEdge),
									   edgeSeq.complement(), edgeId.rc());
		}
	}

	//for each haplotype path, if it's a loop - convert it to a linear edge
	/*GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	for (auto& path : unbranchingPaths)
	{
		if (!path.circular || path.id.strand()) continue;

		for (auto& edge : path.path)
		{
			if (edge->altHaplotype)
			{
				//Logger::get().debug() << "Uncircularized: " << path.edgesStr();
				_graph.disconnectRight(path.path.back());
				_graph.disconnectLeft(_graph.complementEdge(path.path.front()));
				break;
			}
		}
	}*/

	if ((bool)Config::get("remove_alt_edges"))
	{
		std::unordered_set<GraphEdge*> toDelete;
		for (auto& edge: _graph.iterEdges())
		{
			if (edge->altHaplotype) toDelete.insert(edge);
		}
		for (auto& edge : toDelete) _graph.removeEdge(edge);
	}

	_aligner.updateAlignments();
	Logger::get().debug() << "[SIMPL] Collapsed " << numBridged << " haplotypes";
}

DnaSequence HaplotypeResolver::pathSequence(GraphPath& path)
{
	std::string strSeq;
	for (size_t i = 0; i < path.size(); ++i)
	{
		strSeq += _graph.edgeSequences()
			.getSeq(path[i]->seqSegments.front().edgeSeqId).str();
	}
	if (strSeq.empty()) strSeq = "A";
	return DnaSequence(strSeq);
}

void HaplotypeResolver::separeteAdjacentEdges(GraphEdge* inEdge, GraphEdge* outEdge)
{
	GraphNode* newNode = _graph.addNode();

	vecRemove(inEdge->nodeRight->inEdges, inEdge);
	inEdge->nodeRight = newNode;
	newNode->inEdges.push_back(inEdge);

	vecRemove(outEdge->nodeLeft->outEdges, outEdge);
	outEdge->nodeLeft = newNode;
	newNode->outEdges.push_back(outEdge);
}

void HaplotypeResolver::separateDistantEdges(GraphEdge* inEdge, GraphEdge* outEdge,
						  					 EdgeSequence insertSeq, FastaRecord::Id newId)
{
	GraphNode* leftNode = _graph.addNode();
	vecRemove(inEdge->nodeRight->inEdges, inEdge);
	inEdge->nodeRight = leftNode;
	leftNode->inEdges.push_back(inEdge);

	GraphNode* rightNode = _graph.addNode();
	GraphEdge* newEdge = _graph.addEdge(GraphEdge(leftNode, rightNode,
												  newId));
	newEdge->seqSegments.push_back(insertSeq);
	int32_t pathCoverage = (inEdge->meanCoverage +
							outEdge->meanCoverage) / 2;
	newEdge->meanCoverage = pathCoverage;

	vecRemove(outEdge->nodeLeft->outEdges, outEdge);
	outEdge->nodeLeft = rightNode;
	rightNode->outEdges.push_back(outEdge);
}

void HaplotypeResolver::resetEdges()
{
	for (auto& edge : _graph.iterEdges())
	{
		edge->leftLink = nullptr;
		edge->rightLink = nullptr;
		edge->altHaplotype = false;
		edge->altGroupId = -1;
	}
	_bridgingSeqs.clear();
	_nextAltGroupId = 2;
}

//some helper functions for superbubble detection,
//mostly variations of dfs searches
namespace
{
	struct PathWithLen 
	{
		GraphPath path; 
		int length;
	};
	
	//finds any path of length up to maxDepth from the given edge
	//(first and last edges do not count towards length).
	//if there are no paths of that length, return the longest one
	auto anyPath = [](GraphEdge* startEdge, int maxDepth, const RepeatGraph& graph)
	{
		std::vector<PathWithLen> deadEnds;
		std::vector<PathWithLen> queue;
		queue.push_back({{startEdge}, 0});
		
		while (!queue.empty())
		{
			auto curPath = queue.back();
			queue.pop_back();

			//found path that is long enough
			if (curPath.length > maxDepth) return curPath.path;

			bool deadEnd = true;
			for (auto& nextEdge : curPath.path.back()->nodeRight->outEdges)
			{
				bool localRepeat = std::find(curPath.path.rbegin(), curPath.path.rend(),
										  	 nextEdge) != curPath.path.rend();
				bool localComplRepeat = std::find(curPath.path.rbegin(), curPath.path.rend(),
										  	 	  graph.complementEdge(nextEdge)) != curPath.path.rend();
				if (localRepeat || localComplRepeat) continue;
				if (nextEdge->isLooped() && nextEdge->length() < maxDepth) continue;
				deadEnd = false;

				auto newPath = curPath;
				newPath.path.push_back(nextEdge);
				newPath.length += nextEdge->length();
				queue.push_back(newPath);
			}
			if (deadEnd)
			{
				deadEnds.push_back(curPath);
			}

			//heuristic to prevent very slow processing of very tangled components
			const size_t MAX_CANDIDATES = 100;
			if (deadEnds.size() > MAX_CANDIDATES) break;
			if (queue.size() > MAX_CANDIDATES) break;
		}

		//no paths over MAX_DEPTH, return the longest path
		if (deadEnds.empty()) return GraphPath();
		size_t maxId = 0;
		for (size_t i = 1; i < deadEnds.size(); ++i)
		{
			if (deadEnds[i].length > deadEnds[maxId].length) maxId = i;
		}
		return deadEnds[maxId].path;
	};

	//helper structures for dijkstra
	struct DijkstraResult
	{
		DijkstraResult(): failure(false) {}
		std::unordered_map<GraphEdge*, int> dist;
		bool failure;
	};
	struct EdgeWithPriority 
	{
		GraphEdge* edge; 
		int priority;
	};
	auto lowestPri = [](const EdgeWithPriority& e1, const EdgeWithPriority& e2)
		{return e1.priority != e2.priority ? 
				e1.priority > e2.priority : e1.edge > e2.edge;};

	//runs dijskstra algorithm and returns shortest paths
	//lengths from source to all edges. terminates
	//with failure if the distance start to exceed maxBubble
	auto getShortestPathsLen = [](GraphEdge* source, GraphEdge* sink,
								  int maxBubble, DijkstraResult& outRes)
	{
		outRes.dist.clear();
		outRes.failure = false;
		outRes.dist[source] = 0;
		std::set<EdgeWithPriority, decltype(lowestPri)> queue(lowestPri);
		queue.insert({source, 0});

		while (!queue.empty())
		{
			auto curEdge = *queue.begin();
			queue.erase(queue.begin());

			//dead end
			if (curEdge.edge->nodeRight->outEdges.empty()) 
			{
				outRes.failure = true;
				return;
			}

			for (auto& nextEdge : curEdge.edge->nodeRight->outEdges)
			{
				if (nextEdge == sink) continue;

				//we have reached the start, that means there is a bad loop. fail
				if (nextEdge == source)
				{
					outRes.failure = true;
					return;
				}

				int newDist = curEdge.priority + nextEdge->length() + 1;
				if (!outRes.dist.count(nextEdge) || newDist < outRes.dist[nextEdge])
				{
					//path is too long
					if (newDist > maxBubble)
					{
						outRes.failure = true;
						return;
					}

					if (!nextEdge->isLooped())
					{
						queue.erase({nextEdge, outRes.dist[nextEdge]});
						queue.insert({nextEdge, newDist});
					}
					outRes.dist[nextEdge] = newDist;
				}
			}
		}

		outRes.dist.erase(source);
	};

	struct Superbubble
	{
		Superbubble(): start(nullptr), end(nullptr) {}
		GraphEdge* start;
		GraphEdge* end;
		std::unordered_set<GraphEdge*> internalEdges;
		GraphPath refPath;
	};

	Superbubble isRightSuperbubble(GraphEdge* startEdge, int maxBubbleLen,
								   const RepeatGraph& graph, 
								   const std::unordered_set<GraphEdge*> loopedEdges)
	{
		//Logger::get().debug() << "\t\tSearching for ref. path";
		auto refPath = anyPath(startEdge, maxBubbleLen, graph);
		if (refPath.empty()) return Superbubble();

		/*std::string pathStr;
		for (size_t i = 0; i < refPath.size(); ++i)
		{
			pathStr += std::to_string(refPath[i]->edgeId.signedId()) + " -> ";
		}
		Logger::get().debug() << "\t\tReference path: " << pathStr;*/

		for (GraphEdge* endCand : refPath)
		{
			if (endCand == startEdge) continue;
			if (loopedEdges.count(endCand)) continue;
			if (!endCand->nodeLeft->isBifurcation()) continue;
			//if (endCand->nodeLeft->inEdges.size() < 2) continue;

			static DijkstraResult distancesFromSource;
			getShortestPathsLen(startEdge, endCand, maxBubbleLen, 
								distancesFromSource);
			if (distancesFromSource.failure)
			{
				//Logger::get().debug() << "\t\t\tDiverged " << distancesFromSource.dist.size();
				continue;
			}

			/*for (auto& edgeDist : distancesFromSource.dist)
			{
				Logger::get().debug() << "\t\t\tFrom source " << edgeDist.first->edgeId.signedId()
					<< " " << edgeDist.second;
			}*/

			static DijkstraResult distancesFromSink;
			getShortestPathsLen(graph.complementEdge(endCand), 
								graph.complementEdge(startEdge),
								maxBubbleLen, distancesFromSink);
			if (distancesFromSink.failure)
			{
				//Logger::get().debug() << "\t\t\tDiverged " << distancesFromSink.dist.size();
				continue;
			}

			/*for (auto& edgeDist : distancesFromSink.dist)
			{
				Logger::get().debug() << "\t\t\tFrom sink " << edgeDist.first->edgeId.signedId()
					<< " " << edgeDist.second;
			}*/

			//check that set of reachable (internal) edges consistent
			//and each internal edge is connected to source and sink
			//by a path < maxBubble.
			bool goodBubble = true;
			for (auto& edgeDist : distancesFromSource.dist)
			{
				GraphEdge* complEdge = graph.complementEdge(edgeDist.first);
				if (!distancesFromSink.dist.count(complEdge))
				{
					goodBubble = false;
					break;
				}
				int tourLen = edgeDist.second + 
					distancesFromSink.dist[complEdge] - complEdge->length();
				if (tourLen > maxBubbleLen)
				{
					goodBubble = false;
					break;
				}
			}
			for (auto& edgeDist : distancesFromSink.dist)
			{
				GraphEdge* complEdge = graph.complementEdge(edgeDist.first);
				if (!distancesFromSource.dist.count(complEdge))
				{
					goodBubble = false;
					break;
				}
			}

			if (goodBubble)
			{
				Superbubble sb;
				sb.start = startEdge;
				sb.end = endCand;
				sb.refPath = refPath;
				for (auto& edgeDist : distancesFromSource.dist)
				{
					if (edgeDist.first != sb.start && 
						edgeDist.first != sb.end) sb.internalEdges.insert(edgeDist.first);
				}
				//Logger::get().debug() << "Good end: " << sb.end->edgeId.signedId();
				return sb;
			}
		}
		
		return Superbubble();
	}
}

int HaplotypeResolver::findSuperbubbles()
{
	const int MAX_BUBBLE_LEN = Config::get("max_bubble_length");	//50k
	
	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	std::unordered_set<GraphEdge*> loopedEdges;
	for (auto& path : unbranchingPaths)
	{
		if (path.isLooped())
		{
			loopedEdges.insert(path.path.begin(), path.path.end());
		}
	}

	int foundNew = 0;
	std::unordered_set<GraphEdge*> usedEdges;
	for (auto& startEdge : _graph.iterEdges())
	{
		if (loopedEdges.count(startEdge)) continue;
		if (usedEdges.count(startEdge)) continue;
		//if (startEdge->nodeRight->outEdges.size() < 2) continue;
		//if (!startEdge->nodeRight->isBifurcation()) continue;
		
		//require at least two alternative paths (not counting loops)
		//to initiate the bubble. No such requirement for the bubble end though
		int numIn = 0;
		int numOut = 0;
		for (auto& edge : startEdge->nodeRight->inEdges)
		{
			if (!loopedEdges.count(edge)) ++numIn;
		}
		for (auto& edge : startEdge->nodeRight->outEdges)
		{
			if (!loopedEdges.count(edge)) ++numOut;
		}
		//if (numOut < 2) continue;
		//Logger::get().debug() << "\tChecking start: " << startEdge->edgeId.signedId();
		if (numOut < 2 || numIn > 1) continue;

		//if (startEdge->nodeRight->inEdges.size() > 1 ||
		//	startEdge->nodeRight->outEdges.size() < 2) continue;

		//finding superbubble in one direction
		auto fwdBubble = isRightSuperbubble(startEdge, MAX_BUBBLE_LEN,
											_graph, loopedEdges);
		if (!fwdBubble.end || startEdge == fwdBubble.end ||
			startEdge == _graph.complementEdge(fwdBubble.end)) continue;

		//prohibit complex loops to be called as superbubbles
		if (startEdge->nodeRight == fwdBubble.end->nodeLeft) continue;

		//in the opposite direction, both directions must agree
		/*auto revBubble = isRightSuperbubble(_graph.complementEdge(fwdBubble.end),
											MAX_BUBBLE_LEN, _graph, loopedEdges);
		if (!revBubble.end || 
			startEdge != _graph.complementEdge(revBubble.end))
		{
			Logger::get().warning() << "Direction inconsistency! "
				<< startEdge->edgeId.signedId() << " " 
				<< fwdBubble.end->edgeId.signedId();
			continue;
		}*/

		//superbubble found!
		//Logger::get().debug() << "\t\tFound superbubble!";
		usedEdges.insert(_graph.complementEdge(fwdBubble.end));

		//get intermediate edges and mark them
		bool newVariant = true;
		for (auto& edge : fwdBubble.internalEdges)
		{
			if (edge->altHaplotype) newVariant = false;
		}
		if (newVariant) ++foundNew;

		for (auto& edge : fwdBubble.internalEdges)
		{
			edge->altHaplotype = true;
			edge->altGroupId = _nextAltGroupId;
			_graph.complementEdge(edge)->altHaplotype = true;
			_graph.complementEdge(edge)->altGroupId = _nextAltGroupId + 1;
		}
		_nextAltGroupId += 2;

		if (startEdge->rightLink || fwdBubble.end->leftLink) continue;

		//link edges
		_graph.linkEdges(startEdge, fwdBubble.end);
		_graph.linkEdges(_graph.complementEdge(fwdBubble.end), 
						 _graph.complementEdge(startEdge));

		//bridging sequence
		GraphPath bridgePath(fwdBubble.refPath.begin() + 1,
							 std::find(fwdBubble.refPath.begin(), fwdBubble.refPath.end(),
								 	   fwdBubble.end));
		DnaSequence pathSeq = this->pathSequence(bridgePath);
		_bridgingSeqs[std::make_pair(startEdge, fwdBubble.end)] = pathSeq;
		_bridgingSeqs[std::make_pair(_graph.complementEdge(fwdBubble.end), 
									 _graph.complementEdge(startEdge))] = 
												pathSeq.complement();
		//

		Logger::get().debug() << "Superbubble: " << startEdge->edgeDescrLong()
			<< " " << fwdBubble.end->edgeDescrLong();
		std::string pathStr;
		for (auto& edge : fwdBubble.internalEdges)
		{
			pathStr += edge->edgeDescrLong() + " ";
		}
		Logger::get().debug() << "\tInternal: " << pathStr;

		/*if (startEdge->nodeRight->inEdges.size() > 1 ||
			fwdBubble.end->nodeLeft->outEdges.size() > 1)
		if (numIn > 1)
		{
			Logger::get().debug() << " Interesting!";
		}*/
		//
	}

	Logger::get().debug() << "[SIMPL] Masked " << foundNew
		<< " superbubbles";

	return foundNew;
}
