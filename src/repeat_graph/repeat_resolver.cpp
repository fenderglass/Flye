//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cmath>


#include "repeat_resolver.h"
#include "graph_processing.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"

#include <lemon/list_graph.h>
#include <lemon/matching.h>

//Given the path in the graph with a resolved repeat inside,
//separates in into a single unbranching path. The first
//and the last edges of the graphPath parameter
//should correspond to the flanking unique edges
void RepeatResolver::separatePath(const GraphPath& graphPath, 
								  EdgeSequence readSegment, 
								  FastaRecord::Id newId)
{
	//first edge
	GraphNode* leftNode = _graph.addNode();
	vecRemove(graphPath.front()->nodeRight->inEdges, graphPath.front());
	graphPath.front()->nodeRight = leftNode;
	leftNode->inEdges.push_back(graphPath.front());
	int32_t pathCoverage = (graphPath.front()->meanCoverage +
						    graphPath.back()->meanCoverage) / 2;

	//repetitive edges in the middle
	for (size_t i = 1; i < graphPath.size() - 1; ++i)
	{
		graphPath[i]->resolved = true;
		_substractedCoverage[graphPath[i]] += pathCoverage;
		//graphPath[i]->substractedCoverage += pathCoverage;
	}

	GraphNode* rightNode = leftNode;
	if (graphPath.size() > 2)
	{
		rightNode = _graph.addNode();
		GraphEdge* newEdge = _graph.addEdge(GraphEdge(leftNode, rightNode,
													  newId));
		newEdge->seqSegments.push_back(readSegment);
		newEdge->meanCoverage = pathCoverage;
	}

	//last edge
	vecRemove(graphPath.back()->nodeLeft->outEdges, graphPath.back());
	graphPath.back()->nodeLeft = rightNode;
	rightNode->outEdges.push_back(graphPath.back());
}

//Resolves all repeats simulateously through the graph mathcing optimization,
//Given the reads connecting unique edges (or pairs of edges in the transitions graph)
int RepeatResolver::resolveConnections(const std::vector<Connection>& connections, 
									   float minSupport)
{
	ska::flat_hash_map<FastaRecord::Id, 
					   std::vector<const Connection*>> connectIndex;
	for (auto& conn : connections)
	{
		connectIndex[conn.path.front()->edgeId].push_back(&conn);
		connectIndex[conn.path.front()->edgeId.rc()].push_back(&conn);
		connectIndex[conn.path.back()->edgeId].push_back(&conn);
		connectIndex[conn.path.back()->edgeId.rc()].push_back(&conn);
	}

	//Constructs transitions graph using the lemon library
	ska::flat_hash_map<FastaRecord::Id, int> leftCoverage;
	ska::flat_hash_map<FastaRecord::Id, int> rightCoverage;

	ska::flat_hash_map<FastaRecord::Id, int> asmToLemon;
	ska::flat_hash_map<int, FastaRecord::Id> lemonToAsm;
	lemon::ListGraph graph;
	lemon::ListGraph::EdgeMap<int> edgeWeights(graph);

	auto getEdge = [&graph](lemon::ListGraph::Node n1, lemon::ListGraph::Node n2)
	{
		for (lemon::ListGraph::IncEdgeIt edgeIt(graph, n1); 
			 edgeIt != lemon::INVALID; ++edgeIt) 
		{
			if (graph.oppositeNode(n1, edgeIt) == n2) return edgeIt;
		}
		return lemon::ListGraph::IncEdgeIt(lemon::INVALID);
	};

	for (auto& conn : connections)
	{
		GraphEdge* leftEdge = conn.path.front();
		GraphEdge* rightEdge = conn.path.back();

		if (leftEdge->edgeId == rightEdge->edgeId ||
			leftEdge->edgeId == rightEdge->edgeId.rc()) continue;

		++leftCoverage[leftEdge->edgeId];
		++rightCoverage[rightEdge->edgeId.rc()];

		if (!asmToLemon.count(leftEdge->edgeId))
		{
			auto newNode = graph.addNode();
			asmToLemon[leftEdge->edgeId] = graph.id(newNode);
			lemonToAsm[graph.id(newNode)] = leftEdge->edgeId;
		}
		if (!asmToLemon.count(rightEdge->edgeId.rc()))
		{
			auto newNode = graph.addNode();
			asmToLemon[rightEdge->edgeId.rc()] = graph.id(newNode);
			lemonToAsm[graph.id(newNode)] = rightEdge->edgeId.rc();
		}

		auto leftLemonNode = graph.nodeFromId(asmToLemon[leftEdge->edgeId]);
		auto rightLemonNode = graph.nodeFromId(asmToLemon[rightEdge->edgeId.rc()]);
		if (!graph.valid(getEdge(leftLemonNode, rightLemonNode)))
		{
			auto edge = graph.addEdge(leftLemonNode, rightLemonNode);
			edgeWeights[edge] = 0;
		}
		auto edge = getEdge(leftLemonNode, rightLemonNode);
		++edgeWeights[edge];
	}

	//copmutes maximum weight matching on this graph
	lemon::MaxWeightedMatching<lemon::ListGraph> matcher(graph, edgeWeights);
	matcher.run();

	//converting matching to the resolved paths on the graph
	ska::flat_hash_set<FastaRecord::Id> usedEdges;
	std::vector<Connection> uniqueConnections;
	int unresolvedLinks = 0;
	for (auto lemonAsm : lemonToAsm)
	{
		auto mateNode = matcher.mate(graph.nodeFromId(lemonAsm.first));
		if (mateNode == lemon::INVALID) continue;

		FastaRecord::Id leftId = lemonAsm.second;
		FastaRecord::Id rightId = lemonToAsm[graph.id(mateNode)];
		int support = edgeWeights[getEdge(graph.nodeFromId(lemonAsm.first), 
										  mateNode)];

		if (usedEdges.count(leftId)) continue;
		usedEdges.insert(rightId);

		float confidence = (float)support / (leftCoverage[leftId] + 
									  		 rightCoverage[rightId]);

		Logger::get().debug() << "\tConnection " 
			<< leftId.signedId() << "\t" << rightId.rc().signedId()
			<< "\t" << support / 4 << "\t" << confidence;

		if (confidence < minSupport)
		{
			++unresolvedLinks;
			continue;
		}

		std::vector<Connection> spanningConnections;
		for (auto& conn : connectIndex[leftId])
		{
			if ((conn->path.front()->edgeId == leftId && 
				 	conn->path.back()->edgeId == rightId.rc()) ||
				(conn->path.front()->edgeId == rightId && 
				 	conn->path.back()->edgeId == leftId.rc()))
			{
				spanningConnections.push_back(*conn);
			}
		}
		if (spanningConnections.empty())
		{
			Logger::get().warning() << "Empty spanning connections";
			continue;
		}
		std::sort(spanningConnections.begin(), spanningConnections.end(),
				  [](const Connection c1, const Connection c2)
				{return c1.readSeq.length() < c2.readSeq.length();});
		uniqueConnections
			.push_back(spanningConnections[spanningConnections.size() / 2]);
	}

	//separates the resolved paths in the graph
	for (auto& conn : uniqueConnections)
	{
		FastaRecord::Id edgeId = _graph.newEdgeId();

		std::string description = "edge_" + std::to_string(edgeId.signedId()) + 
				"_0_" + _readSeqs.getRecord(conn.readSeq.readId).description + "_" +
				std::to_string(conn.readSeq.start) + "_" + 
				std::to_string(conn.readSeq.end);
		EdgeSequence edgeSeq = 
			_graph.addEdgeSequence(_readSeqs.getSeq(conn.readSeq.readId),
								   conn.readSeq.start, conn.readSeq.length(),
								   description);

		this->separatePath(conn.path, edgeSeq, edgeId);
		this->separatePath(_graph.complementPath(conn.path), 
						   edgeSeq.complement(), edgeId.rc());
	}

	Logger::get().debug() << "Resolved: " << uniqueConnections.size() << " links: "
						  << connections.size() / 2;
	Logger::get().debug() << "Unresolved: " << unresolvedLinks;

	return uniqueConnections.size();
}

bool RepeatResolver::checkByReadExtension(const GraphEdge* checkEdge,
										  const std::vector<GraphAlignment>& alignments)
{
	ska::flat_hash_map<GraphEdge*, std::vector<int>> outFlanks;
	ska::flat_hash_map<GraphEdge*, std::vector<int>> outSpans;
	int lowerBound = 0;
	for (auto& aln : alignments)
	{ 
		bool passedStart = false;
		int leftFlank = 0;
		int leftCoord = 0;
		bool foundUnique = false;
		for (size_t i = 0; i < aln.size(); ++i)
		{
			 //only high quality flanking alignments
			//if (aln[i].overlap.seqDivergence > 0.15) continue;

			if (!passedStart && aln[i].edge == checkEdge)
			{
				passedStart = true;
				leftFlank = aln[i].overlap.curEnd - aln[0].overlap.curBegin;
				leftCoord = aln[i].overlap.curEnd;
				continue;
			}
			if (passedStart && !aln[i].edge->repetitive)
			{
				if (aln[i].edge->edgeId != checkEdge->edgeId &&
					aln[i].edge->edgeId != checkEdge->edgeId.rc())
				{
					int rightFlank = aln.back().overlap.curEnd -
									 aln[i].overlap.curBegin;
					int alnSpan = aln[i].overlap.curBegin - leftCoord;
					outFlanks[aln[i].edge].push_back(std::min(leftFlank, rightFlank));
					outSpans[aln[i].edge].push_back(alnSpan);
				}
				foundUnique = true;
				break;
			}
		}
		if (!foundUnique)
		{
			lowerBound = std::max(lowerBound, 
								  aln.back().overlap.curBegin - leftCoord);
		}
	}

	//check if there is agreement
	int maxSupport = 0;
	for (auto& outConn : outFlanks)
	{
		if (maxSupport < (int)outConn.second.size())
		{
			maxSupport = outConn.second.size();
		}
	}

	int uniqueMult = 0;
	int minSupport = maxSupport / (int)Config::get("out_paths_ratio");
	//if there is at least one extension supported by more than 1 read,
	//make minimum support at least 1
	if (maxSupport > 1) minSupport = std::max(minSupport, 1);

	for (auto& outConn : outFlanks) 
	{
		if ((int)outConn.second.size() > minSupport)
		{
			++uniqueMult;
		}
	}
	
	if (uniqueMult > 1) 
	{
		Logger::get().debug() << "Starting " 
			<< checkEdge->edgeId.signedId() << " aln:" << alignments.size()
			<< " minSpan:" << lowerBound;
		for (auto& outEdgeCount : outFlanks)
		{
			int maxFlank = *std::max_element(outEdgeCount.second.begin(),
											 outEdgeCount.second.end());
			int minSpan = *std::min_element(outSpans[outEdgeCount.first].begin(),
											outSpans[outEdgeCount.first].end());

			std::string star = outEdgeCount.first->repetitive ? "R" : " ";
			std::string loop = outEdgeCount.first->isLooped() ? "L" : " ";
			std::string tip = outEdgeCount.first->isTip() ? "T" : " ";
			Logger::get().debug() << "\t" << star << " " << loop << " " << tip << " "
				<< outEdgeCount.first->edgeId.signedId() << "\tnum:" << outEdgeCount.second.size()
				<< "\tflank:" << maxFlank << "\tspan:" << minSpan;
		}
		return true;
	}
	return false;
}

//Classifies all edges into unique and repetitive based on the coverage + 
//alignment information - one of the key steps here.
void RepeatResolver::findRepeats()
{
	Logger::get().debug() << "Finding repeats";

	ska::flat_hash_map<GraphEdge*, 
					   std::vector<GraphAlignment>> alnIndex;
	for (auto& aln : _aligner.getAlignments())
	{
		if (aln.size() > 1)
		{
			for (auto& edgeAln : aln)
			{
				alnIndex[edgeAln.edge].push_back(aln);
			}
		}
	}

	//all edges are unique at the beginning
	for (auto& edge : _graph.iterEdges())
	{
		edge->repetitive = false;
	}

	//Will operate on unbranching paths rather than single edges
	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	ska::flat_hash_map<FastaRecord::Id, UnbranchingPath*> idToPath;
	for (auto& path : unbranchingPaths) idToPath[path.id] = &path;
	auto complPath = [&idToPath](UnbranchingPath* path)
	{
		if (idToPath.count(path->id.rc()))
		{
			return idToPath[path->id.rc()];
		}
		return path;	//self-complement
	};
	auto markRepetitive = [](UnbranchingPath* path)
	{
		for (auto& edge : path->path) edge->repetitive = true;
	};

	for (auto& path : unbranchingPaths)
	{
		if (!path.id.strand()) continue;

		//mark paths with high coverage as repetitive
		if (!Parameters::get().unevenCoverage)
		{
			if (path.meanCoverage > _multInf.getUniqueCovThreshold())
			{
				markRepetitive(&path);
				markRepetitive(complPath(&path));
				Logger::get().debug() << "Cov: " 
					<< path.edgesStr() << "\t" << path.length << "\t" 
					<< path.meanCoverage;
			}
		}

		//self-complements are repetitive
		for (auto& edge : path.path)
		{
			if (edge->selfComplement)
			{
				Logger::get().debug() << "Self-compl: " << path.edgesStr();
				markRepetitive(&path);
				markRepetitive(complPath(&path));
				break;
			}
		}

		//This was kinda hacky originally - probably should be removed now
		//with all other repeat detection improvements
		//tandem repeats
		/*if (path.path.size() == 1 && path.path.front()->isLooped())
		{
			ska::flat_hash_set<FastaRecord::Id> seen;
			for (auto& seg : path.path.front()->seqSegments)
			{
				if (seen.count(seg.origSeqId))
				{
					markRepetitive(&path);
					markRepetitive(complPath(&path));
					Logger::get().debug() << "Tan: " 
						<< path.edgesStr() << "\t" << path.length << "\t" 
						<< path.meanCoverage;
					break;
				}
				seen.insert(seg.origSeqId);
			}
		}*/
	}

	//Finally, using the read alignments
	//order might be important, process short edges first
	std::vector<UnbranchingPath*> sortedPaths;
	for (auto& path : unbranchingPaths) sortedPaths.push_back(&path);
	std::sort(sortedPaths.begin(), sortedPaths.end(),
			  [](const UnbranchingPath* p1, const UnbranchingPath* p2) 
			  {return p1->length < p2->length;});

	//in the case of metagenome do 2 passes, since some small
	//edges might not be detected from the 1st iteration
	//if tey are partes of mosaic repeats. In the case of
	//uniform coverage, this edges are typically detected using coverage
	size_t numIters = !Parameters::get().unevenCoverage ? 1 : 2;
	for (size_t i = 0; i < numIters; ++i)
	{
		Logger::get().debug() << "Repeat detection iteration " << i + 1;
		for (auto& path : sortedPaths)
		{
			if (!path->id.strand()) continue;
			if (path->path.front()->repetitive) continue;
			//if (path->length > (int)Config::get("unique_edge_length")) continue;

			bool rightRepeat = 
				this->checkByReadExtension(path->path.back(), 
										   alnIndex[path->path.back()]);
			bool leftRepeat = 
				this->checkByReadExtension(complPath(path)->path.back(), 
										   alnIndex[complPath(path)->path.back()]);
			if (rightRepeat || leftRepeat)
			{
				markRepetitive(path);
				markRepetitive(complPath(path));
				
				Logger::get().debug() << "Mult: " 
					<< path->edgesStr() << "\t" << path->length << "\t" 
					<< path->meanCoverage << "\t" " ("
					<< leftRepeat << "," << rightRepeat << ")";
			}
		}
	}

	//now, check for this structure >-<, in case read alignments were not enough
	/*for (auto& path : sortedPaths)
	{
		if (path->path.front()->repetitive || path->isLooped()) continue;
		if (path->path.front()->nodeLeft->outEdges.size() > 1 ||
			path->path.back()->nodeRight->inEdges.size() > 1) continue;

		int numIn = 0;
		int numOut = 0;
		for (auto& edge: path->path.front()->nodeLeft->inEdges)
		{
			if (!edge->repetitive) ++numIn;
		}
		for (auto& edge: path->path.back()->nodeRight->outEdges)
		{
			if (!edge->repetitive) ++numOut;
		}
		if (numIn > 1 && numOut > 1)
		{
			Logger::get().debug() << "Structure: " << path->edgesStr();
			markRepetitive(path);
			markRepetitive(complPath(path));
		}
	}*/
}

void RepeatResolver::finalizeGraph()
{
	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	for (auto& path : unbranchingPaths)
	{
		if (!path.id.strand()) continue;

		bool highCoverage = (float)path.meanCoverage > 
							_multInf.getUniqueCovThreshold();

		if (!path.path.front()->selfComplement &&
			path.path.front()->repetitive &&
			path.length > (int)Config::get("unique_edge_length") &&
			(Parameters::get().unevenCoverage || !highCoverage))
		{
			for (auto& edge : path.path)
			{
				edge->repetitive = false;
				_graph.complementEdge(edge)->repetitive = false;
			}

			Logger::get().debug() << "Fixed: " 
				<< path.edgesStr() << "\t" << path.length << "\t" 
				<< path.meanCoverage;
		}
	}

	//apply coverage substractions that were made during repeat resolution
	for (auto& path : unbranchingPaths)
	{
		if (path.isLooped()) continue;
		for (auto& edge : path.path)
		{
			edge->meanCoverage = std::max(0, (int)edge->meanCoverage - 
											  _substractedCoverage[edge]);
		}
	}
}

//Iterates repeat detection and resolution until
//no new repeats are resolved
void RepeatResolver::resolveRepeats()
{
	//make the first iteration resolve only high-confidence connections
	while (true)
	{
		float minSupport = Config::get("min_repeat_res_support");
		auto connections = this->getConnections();
		int resolvedConnections = 
			this->resolveConnections(connections, minSupport);

		this->clearResolvedRepeats();
		
		if (!resolvedConnections) break;
		this->findRepeats();
	}

	GraphProcessor proc(_graph, _asmSeqs);
	proc.fixChimericJunctions();
	_aligner.updateAlignments();
}


//extracts conenctions between pairs of unique edges from
//read alignments
std::vector<RepeatResolver::Connection> 
	RepeatResolver::getConnections()
{
	
	auto safeEdge = [this](GraphEdge* edge)
	{
		return !edge->isRepetitive();
	};

	int totalSafe = 0;
	for (GraphEdge* edge : _graph.iterEdges())
	{
		if (edge->edgeId.strand() && safeEdge(edge)) ++totalSafe;
	}
	Logger::get().debug() << "Total unique edges: " << totalSafe;

	const int32_t MAGIC_100 = 100;
	std::vector<Connection> readConnections;
	for (auto& readPath : _aligner.getAlignments())
	{
		GraphAlignment currentAln;
		int32_t readStart = 0;
		for (auto& aln : readPath)
		{
			if (currentAln.empty()) 
			{
				if (!safeEdge(aln.edge)) continue;
				readStart = aln.overlap.curEnd + aln.overlap.extLen - 
							aln.overlap.extEnd;
				readStart = std::min(readStart, aln.overlap.curLen - MAGIC_100);
			}

			currentAln.push_back(aln);
			if (safeEdge(aln.edge) && currentAln.front().edge != aln.edge)
			{
				if (!currentAln.back().edge->nodeLeft->isBifurcation() &&
					!currentAln.front().edge->nodeRight->isBifurcation()) continue;

				//if (currentAln.front().overlap.seqDivergence > 0.15 ||
				//	currentAln.back().overlap.seqDivergence > 0.15) continue;
				
				int32_t flankScore = std::min(currentAln.front().overlap.curRange(),
											  currentAln.back().overlap.curRange());
				GraphPath currentPath;
				for (auto& aln : currentAln) currentPath.push_back(aln.edge);
				GraphPath complPath = _graph.complementPath(currentPath);

				int32_t readEnd = aln.overlap.curBegin - aln.overlap.extBegin;

				//TODO: fix this ad-hoc fix. Currently, if read connects
				//two consecutive edges (for example, when resolving chimera junctions,
				//we still would insert a tiny bit of read sequence as a placeholder.
				//Probably, wouldn't harm, but who knows..
				readEnd = std::max(readStart + MAGIC_100 - 1, readEnd);	
				if (readStart < 0 || readEnd >= aln.overlap.curLen)
				{
					Logger::get().warning() 
						<< "Something is wrong with bridging read sequence";
					//Logger::get().warning() << readStart << " " 
					//	<< readEnd << " " << aln.overlap.curLen;
					break;
				}

				/*std::string description = "read_seq_" + 
					std::to_string(aln.overlap.curId.signedId());
				EdgeSequence edgeSeq = 
					_graph.addEdgeSequence(_readSeqs.getSeq(aln.overlap.curId),
										   readStart, readEnd - readStart,
										   description);*/
				ReadSequence readSeq = {aln.overlap.curId, readStart, readEnd};
				ReadSequence complRead = {aln.overlap.curId.rc(), 
										  aln.overlap.curLen - readEnd - 1,
										  aln.overlap.curLen - readStart - 1};
				readConnections.push_back({currentPath, readSeq, flankScore});
				readConnections.push_back({complPath, complRead, flankScore});

				currentAln.clear();
				currentAln.push_back(aln);
				readStart = aln.overlap.curEnd + aln.overlap.extLen - 
							aln.overlap.extEnd;
				readStart = std::min(readStart, aln.overlap.curLen - MAGIC_100);
			}
		}
	}

	return readConnections;
}

//cleans up the graph after repeat resolution
void RepeatResolver::clearResolvedRepeats()
{
	//const int MIN_LOOP = Parameters::get().minimumOverlap;
	auto nextEdge = [](GraphNode* node)
	{
		for (auto edge : node->outEdges)
		{
			if (!edge->isLooped()) return edge;
		}
		return (GraphEdge*)nullptr;
	};

	auto shouldRemove = [](GraphEdge* edge)
	{
		//return edge->isRepetitive() && edge->resolved;
		return edge->resolved;
	};

	ska::flat_hash_set<GraphNode*> toRemove;

	for (auto& node : _graph.iterNodes())
	{
		//separated nodes
		if (node->neighbors().size() == 0)
		{
			bool resolved = true;
			for (auto& edge : node->outEdges) 
			{
				if (!shouldRemove(edge)) resolved = false;
			}

			if (resolved) toRemove.insert(node);
		}

		//other nodes
		if (!node->isEnd()) continue;

		GraphEdge* direction = nextEdge(node);
		if (!direction) continue;

		GraphPath traversed;
		traversed.push_back(direction);
		GraphNode* curNode = direction->nodeRight;
		while (curNode->isResolved())
		{
			traversed.push_back(nextEdge(curNode));
			curNode = traversed.back()->nodeRight;
		}
		if (traversed.empty()) continue;

		bool removeLast = curNode->isEnd();
		bool resolvedRepeat = true;
		for (auto& edge : traversed) 
		{
			if (!shouldRemove(edge)) resolvedRepeat = false;
		}

		GraphPath complPath = _graph.complementPath(traversed);
		if (resolvedRepeat)
		{
			//first-last
			toRemove.insert(traversed.front()->nodeLeft);
			if (removeLast) toRemove.insert(complPath.front()->nodeLeft);

			//middle nodes
			for (size_t i = 0; i < traversed.size() - 1; ++i)
			{
				toRemove.insert(traversed[i]->nodeRight);
				toRemove.insert(complPath[i]->nodeRight);
			}

			//last-first
			if (removeLast) toRemove.insert(traversed.back()->nodeRight);
			toRemove.insert(complPath.back()->nodeRight);
		}
	}

	for (auto node : toRemove) _graph.removeNode(node);
	_aligner.updateAlignments();
}
