//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "multiplicity_inferer.h"
#include "graph_processing.h"
#include "../common/disjoint_set.h"
#include "../common/utils.h"
#include <cmath>



//Estimates the mean coverage and assingns edges multiplicity accordingly
void MultiplicityInferer::estimateCoverage()
{
	const int WINDOW = Config::get("coverage_estimate_window");

	//alternative coverage
	std::unordered_map<GraphEdge*, std::vector<int32_t>> wndCoverage;

	for (auto& edge : _graph.iterEdges())
	{
		size_t numWindows = edge->length() / WINDOW;
		wndCoverage[edge].assign(numWindows, 0);
	}

	for (auto& path : _aligner.getAlignments())
	{
		for (size_t pathId = 0; pathId < path.size(); ++pathId)
		{
			auto& edgeCov = wndCoverage[path[pathId].edge];
			int covFrom = std::max(0, path[pathId].overlap.extBegin / WINDOW + 1);
			int covTo = std::min((int)edgeCov.size(), path[pathId].overlap.extEnd / WINDOW);

			//for intermediae alignments, cover the entire edge
			if (pathId > 0) covFrom = 0;
			if (pathId < path.size() - 1) covTo = edgeCov.size();

			for (int i = covFrom; i < covTo; ++i) ++edgeCov[i];
		}
	}

	int64_t sumCov = 0;
	int64_t sumLength = 0;
	for (auto& edgeCoverage : wndCoverage)
	{
		for (auto& cov : edgeCoverage.second)
		{
			sumCov += (int64_t)cov;
			++sumLength;
		}
	}
	_meanCoverage = (sumLength != 0) ? sumCov / sumLength : /*defaut*/ 1;

	Logger::get().info() << "Mean edge coverage: " << _meanCoverage;

	std::vector<int32_t> edgesCoverage;
	for (auto edge : _graph.iterEdges())
	{
		if (wndCoverage[edge].empty()) continue;

		GraphEdge* complEdge = _graph.complementEdge(edge);
		int32_t medianCov = (median(wndCoverage[edge]) + 
						 	 median(wndCoverage[complEdge])) / 2;

		int estMult = std::round((float)medianCov / _meanCoverage);
		if (estMult == 1)
		{
			edgesCoverage.push_back(medianCov);
		}

		//std::string match = estMult != edge->multiplicity ? "*" : " ";
		std::string covStr;

		Logger::get().debug() << edge->edgeId.signedId() << "\tlen:"
				<< edge->length() << "\tcov:" << medianCov << "\tmult:"
				<< (float)medianCov / _meanCoverage;

		//edge->multiplicity = estMult;
		edge->meanCoverage = medianCov;
	}

	_uniqueCovThreshold = /*default*/ 2;
	if (!edgesCoverage.empty())
	{
		const float MULT = (float)Config::get("repeat_edge_cov_mult");	//1.75
		_uniqueCovThreshold = MULT * quantile(edgesCoverage, 75);
	}
	Logger::get().debug() << "Unique coverage threshold " << _uniqueCovThreshold;
}

int MultiplicityInferer::resolveForks()
{
	//const int UNIQUE_LEN = (int)Config::get("unique_edge_length");
	const int MAJOR_TO_MINOR = (int)Config::get("weak_detach_rate");

	int numDisconnected = 0;
	std::vector<GraphNode*> originalNodes(_graph.iterNodes().begin(), 
									      _graph.iterNodes().end());
	for (GraphNode* node : originalNodes)
	{
		//only forks: one in, two out
		if (node->inEdges.size() != 1 ||
			node->outEdges.size() != 2) continue;

		GraphEdge* inEdge = node->inEdges.front();
		GraphEdge* outMajor = node->outEdges[0];
		GraphEdge* outMinor = node->outEdges[1];
		if (outMinor->meanCoverage > outMajor->meanCoverage)
		{
			std::swap(outMajor, outMinor);
		}

		//confirming the correct structure
		if (inEdge->selfComplement || inEdge->isLooped() ||
			outMajor->selfComplement || outMajor->isLooped() ||
			outMinor->selfComplement || outMinor->isLooped()) continue;

		//we want out input edge to be unique. This is not the
		//most reliable way to tell, but at least something
		//if (inEdge->length() < UNIQUE_LEN) continue;

		//we want coverage of major edges significantly higher than minor
		if (std::min(outMajor->meanCoverage, inEdge->meanCoverage) < 
			outMinor->meanCoverage * MAJOR_TO_MINOR) continue;

		//looks like all is good
		Logger::get().debug() << "Disconnected fork: " << outMinor->edgeId.signedId();
		_graph.disconnectLeft(outMinor);
		_graph.disconnectRight(_graph.complementEdge(outMinor));
		++numDisconnected;
	}

	_aligner.updateAlignments();
	Logger::get().debug() << "[SIMPL] Simplified " << numDisconnected << " forks";
	return numDisconnected;
}


//removes the edges with low coverage, with an option to
//only remove tips
int MultiplicityInferer::removeUnsupportedEdges(bool onlyTips)
{
	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	int32_t coverageThreshold = 0;
	const int MIN_CUTOFF = std::round((float)Config::get("min_read_cov_cutoff"));
	if (!Parameters::get().unevenCoverage)
	{
		coverageThreshold = std::round((float)this->getMeanCoverage() / 
										Config::get("graph_cov_drop_rate"));
		coverageThreshold = std::max(MIN_CUTOFF, coverageThreshold);
	}
	else
	{
		coverageThreshold = MIN_CUTOFF;
	}
	Logger::get().debug() << "Read coverage cutoff: " << coverageThreshold;

	std::unordered_set<GraphEdge*> toRemove;
	int removedPaths = 0;
	for (auto& path : unbranchingPaths)
	{
		if (!path.id.strand()) continue;

		//check if it's a tip
		if (onlyTips && !path.path.back()->isRightTerminal()) continue;

		if (path.meanCoverage < coverageThreshold)
		{
			++removedPaths;
			for (auto& edge : path.path)
			{
				toRemove.insert(edge);
				toRemove.insert(_graph.complementEdge(edge));
			}
		}
	}

	for (auto& edge : toRemove) _graph.removeEdge(edge);
	Logger::get().debug() << "[SIMPL] Removed " << removedPaths
		<< " paths with low coverage";

	_aligner.updateAlignments();
	return toRemove.size() / 2;
}

int MultiplicityInferer::disconnectMinorPaths()
{
	const int DETACH_RATE = (int)Config::get("weak_detach_rate");
	const int MAX_LEN = 50000;

	auto nodeDegree = [](GraphNode* node)
	{
		std::vector<int> coverages;
		for (auto& edge : node->inEdges) 
		{
			if (!edge->isLooped()) coverages.push_back(edge->meanCoverage);
		}
		for (auto& edge : node->outEdges) 
		{
			if (!edge->isLooped()) coverages.push_back(edge->meanCoverage);
		}
		if (coverages.size() < 3) return 0;
		return median(coverages);

		/*int maxIn = 0;
		int maxOut = 0;
		for (auto& edge : node->inEdges) 
		{
			if (!edge->isLooped()) maxIn = std::max(maxIn, edge->meanCoverage);
		}
		for (auto& edge : node->outEdges) 
		{
			if (!edge->isLooped()) maxOut = std::max(maxOut, edge->meanCoverage);
		}
		return std::min(maxIn, maxOut);*/
	};

	int numDisconnected = 0;
	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	std::unordered_set<FastaRecord::Id> toRemove;
	for (auto& path : unbranchingPaths)
	{
		if (!path.id.strand() || 
			path.isLooped() ||
			path.path.front()->selfComplement ||
			path.length > MAX_LEN) continue;
		if (path.nodeLeft()->inEdges.empty() ||
			path.nodeRight()->outEdges.empty())	continue; //already detached or tip

		bool weakLeft = path.nodeLeft()->inEdges.empty() || 
						nodeDegree(path.nodeLeft()) > path.meanCoverage * DETACH_RATE;
		bool weakRight = path.nodeRight()->outEdges.empty() || 
						 nodeDegree(path.nodeRight()) > path.meanCoverage * DETACH_RATE;
		if (weakLeft && weakRight) toRemove.insert(path.id);
	}
	
	for (auto& path : unbranchingPaths)
	{
		if (toRemove.count(path.id))
		{
			_graph.disconnectLeft(path.path.front());
			_graph.disconnectLeft(_graph.complementEdge(path.path.back()));
			_graph.disconnectRight(path.path.back());
			_graph.disconnectRight(_graph.complementEdge(path.path.front()));
			++numDisconnected;
			Logger::get().debug() << "Fragile path: " << path.edgesStr();
		}
	}

	_aligner.updateAlignments();
	Logger::get().debug() << "[SIMPL] Disconnected "
		<< numDisconnected << " minor paths";

	return numDisconnected;
}

//Checks each node in the graph if all edges form a single
//connectivity cluster based on read alignment. If
//there are multiple cluster, the node is split into 
//multiple corresponding nodes. This, for example,
//addresses chimeric connections
int MultiplicityInferer::splitNodes()
{
	static const int MIN_JCT_SUPPORT = 1;

	Logger::get().debug() << "Splitting nodes";
	int numSplit = 0;

	//storing connectivity information
	std::unordered_map<GraphEdge*, 
					   std::unordered_map<GraphEdge*, int>> readSupport;
	for (auto& readPath : _aligner.getAlignments())
	{
		if (readPath.size() < 2) continue;
		
		for (size_t i = 0; i < readPath.size() - 1; ++i)
		{
			//if (readPath[i].edge == readPath[i + 1].edge &&
			//	readPath[i].edge->isLooped()) continue;
			if (readPath[i].edge->edgeId == 
				readPath[i + 1].edge->edgeId.rc()) continue;

			++readSupport[readPath[i].edge][readPath[i + 1].edge];
		}
	}

	std::unordered_set<GraphNode*> usedNodes;
	std::vector<GraphNode*> originalNodes(_graph.iterNodes().begin(), 
									      _graph.iterNodes().end());
	for (auto& nodeToSplit : originalNodes)
	{
		if (nodeToSplit->inEdges.size() < 2 ||
			nodeToSplit->outEdges.size() < 2) continue;
		if (usedNodes.count(nodeToSplit)) continue;
		usedNodes.insert(_graph.complementNode(nodeToSplit));
		bool selfComplNode = nodeToSplit == _graph.complementNode(nodeToSplit);

		//initializing sets (to cluster them later)
		struct EdgeDir
		{
			GraphEdge* edge;
			bool isInput;
		};
		typedef SetNode<EdgeDir> SetElement;
		SetVec<EdgeDir> allElements;
		std::unordered_map<GraphEdge*, SetElement*> inputElements;
		std::unordered_map<GraphEdge*, SetElement*> outputElements;
		for (GraphEdge* edge : nodeToSplit->inEdges) 
		{
			allElements.push_back(new SetElement({edge, true}));
			inputElements[edge] = allElements.back();
		}
		for (GraphEdge* edge : nodeToSplit->outEdges) 
		{
			allElements.push_back(new SetElement({edge, false}));
			outputElements[edge] = allElements.back();
		}

		//grouping edges if they are connected by reads
		for (GraphEdge* inEdge : nodeToSplit->inEdges)
		{
			for (auto outEdge : readSupport[inEdge])
			{
				if (outEdge.second >= MIN_JCT_SUPPORT)
				{
					unionSet(inputElements[inEdge], 
							 outputElements[outEdge.first]);
				}
			}
		}

		auto clusters = groupBySet(allElements);
		if (clusters.size() > 1)	//need to split the node!
		{
			numSplit += 1;
			//Logger::get().debug() << "Node " 
			//	<< nodeToSplit->inEdges.size() + nodeToSplit->outEdges.size()
			//	<< " clusters: " << clusters.size() << " " << selfComplNode;

			//for (auto& cl : clusters)
			//{
				//Logger::get().debug() << "\tCl: " << cl.second.size();
				//for (auto edgeDir : cl.second)
				//{
				//	Logger::get().debug() << "\t\t" << edgeDir.edge->edgeId.signedId() << " " 
				//		<< edgeDir.edge->length() << " " << edgeDir.edge->meanCoverage << " "
				//		<< edgeDir.isInput;
				//}
			//}

			for (auto& cl : clusters)
			{
				auto switchNode = [](GraphEdge* edge, 
									 GraphNode* newNode,
									 bool isInput)
				{
					if (!isInput)
					{
						vecRemove(edge->nodeLeft->outEdges, edge);
						edge->nodeLeft = newNode;
						newNode->outEdges.push_back(edge);
					}
					else
					{
						vecRemove(edge->nodeRight->inEdges, edge);
						edge->nodeRight = newNode;
						newNode->inEdges.push_back(edge);
					}

				};

				GraphNode* newNode = _graph.addNode();
				GraphNode* newComplNode = _graph.addNode();
				for (auto edgeDir : cl.second)
				{
					GraphEdge* complEdge = _graph.complementEdge(edgeDir.edge);
					//GraphNode* complSplit = _graph.complementNode(nodeToSplit);
					switchNode(edgeDir.edge, newNode, edgeDir.isInput);
					//if (!edgeDir.edge->selfComplement && !selfComplNode)
					if (!selfComplNode)
					{
						switchNode(complEdge, newComplNode, 
								   !edgeDir.isInput);
					}
				}
			}
		}
	}

	_aligner.updateAlignments();
	Logger::get().debug() << "[SIMPL] Split " << numSplit << " nodes";

	return numSplit;
}


//Disconnects edges, which had low number of reads that connect them
//with the rest of the graph. #of reads is relative to the
//edge coverage
int MultiplicityInferer::removeUnsupportedConnections()
{
	//static const int MIN_JCT_SUPPORT = 2;

	std::unordered_map<GraphEdge*, int32_t> rightConnections;
	std::unordered_map<GraphEdge*, int32_t> leftConnections;

	for (auto& readPath : _aligner.getAlignments())
	{
		if (readPath.size() < 2) continue;

		for (size_t i = 0; i < readPath.size() - 1; ++i)
		{
			//if (readPath[i].edge == readPath[i + 1].edge &&
			//	readPath[i].edge->isLooped()) continue;
			if (readPath[i].edge->edgeId == 
				readPath[i + 1].edge->edgeId.rc()) continue;

			++rightConnections[readPath[i].edge];
			++leftConnections[readPath[i + 1].edge];
			GraphEdge* complLeft = _graph.complementEdge(readPath[i].edge);
			GraphEdge* complRight = _graph.complementEdge(readPath[i + 1].edge);
			++rightConnections[complRight];
			++leftConnections[complLeft];
		}
	}

	int numDisconnected = 0;
	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand() || edge->isLooped()) continue;
		GraphEdge* complEdge = _graph.complementEdge(edge);

		//int32_t coverageThreshold = edge->meanCoverage / 
		//						Config::get("graph_cov_drop_rate");
		//coverageThreshold = std::max(MIN_JCT_SUPPORT, coverageThreshold);
		//int32_t coverageThreshold = MIN_JCT_SUPPORT;
		int32_t coverageThreshold = (edge->meanCoverage >= 10) ? 2 : 1;

		//Logger::get().debug() << "Adjacencies: " << edge->edgeId.signedId() << " "
		//	<< leftConnections[edge] / 2 << " " << rightConnections[edge] / 2;

		if (!edge->nodeRight->isEnd() &&
			edge->nodeRight->isBifurcation() &&
			rightConnections[edge] / 2 < coverageThreshold)
		{
			++numDisconnected;
			Logger::get().debug() << "Chimeric right: " <<
				edge->edgeId.signedId() << " " << rightConnections[edge] / 2;

			_graph.disconnectRight(edge);
			_graph.disconnectLeft(complEdge);

			if (edge->selfComplement) continue;	//already discinnected
		}
		if (!edge->nodeLeft->isEnd() &&
			edge->nodeLeft->isBifurcation() &&
			leftConnections[edge] / 2 < coverageThreshold)
		{
			++numDisconnected;
			Logger::get().debug() << "Chimeric left: " <<
				edge->edgeId.signedId() << " " << leftConnections[edge] / 2;

			_graph.disconnectLeft(edge);
			_graph.disconnectRight(complEdge);
		}
	}

	Logger::get().debug() << "[SIMPL] Disconnected " << numDisconnected << " edges";

	_aligner.updateAlignments();
	return numDisconnected;
}

void MultiplicityInferer::trimTipsIteration(int& outShort, int& outLong)
{
	const int SHORT_TIP = Config::get("short_tip_length");
	const int LONG_TIP = Config::get("long_tip_length");
	const int COV_RATE = (int)Config::get("tip_coverage_rate");
	const int LEN_RATE = (int)Config::get("tip_length_rate");

	std::unordered_set<FastaRecord::Id> toRemove;
	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	std::unordered_map<GraphEdge*, UnbranchingPath*> ubIndex;
	for (auto& path : unbranchingPaths)
	{
		ubIndex[path.path.front()] = &path;
		ubIndex[path.path.back()] = &path;
		//for (auto& edge: path.path) ubIndex[edge] = &path;
	}

	int shortClipped = 0;
	int longClipped = 0;

	for (auto& tipPath : unbranchingPaths)
	{
		if (!tipPath.path.back()->isRightTerminal()) continue;	//right-tip
		if (tipPath.nodeLeft()->outEdges.size() == 1) continue;	//already detached from the left
		if (tipPath.path.front()->selfComplement) continue;		//never-ever!

		//short tip, remove regardless of coverage
		if (tipPath.length < SHORT_TIP)
		{
			toRemove.insert(tipPath.id);
			++shortClipped;
			continue;
		}

		//longer then max long tip, continue
		if (tipPath.length > LONG_TIP) continue;

		//tip longer than short, and shorter than long :)
		//need to check the graph structure and coverage
		
		//get "true path" entrance and exit. There must be
		//exactly one of each. True edges should be of sufficient
		//length and/or continue into the rest of the graph
		GraphNode* tipNode = tipPath.nodeLeft();
		std::vector<UnbranchingPath*> entrances;
		for (GraphEdge* edge : tipNode->inEdges)
		{
			UnbranchingPath& path = *ubIndex[edge];
			if (path.path.back() == edge)
			{
				if (path.length > LEN_RATE * tipPath.length ||
					path.nodeLeft()->inEdges.size() > 0) entrances.push_back(&path);
			}
		}
		std::vector<UnbranchingPath*> exits;
		for (GraphEdge* edge : tipNode->outEdges)
		{
			UnbranchingPath& path = *ubIndex[edge];
			if (path.path.front() == edge)
			{
				if (path.length > LEN_RATE * tipPath.length ||
					path.nodeRight()->outEdges.size() > 0) exits.push_back(&path);
			}
		}
		if (entrances.size() != 1 || exits.size() != 1) continue;

		//remove the tip if its coverage or length is
		//significantly lower than the true path's
		int trueCov = std::max(entrances.front()->meanCoverage, 
					 		   exits.front()->meanCoverage);
		int trueLen = std::max(entrances.front()->length, exits.front()->length);
		if (trueCov > COV_RATE * tipPath.meanCoverage ||
			trueLen > LEN_RATE * tipPath.length)
		{
			toRemove.insert(tipPath.id);
			++longClipped;
		}
	}
	
	for (auto& path : unbranchingPaths)
	{
		if (toRemove.count(path.id))
		{
			Logger::get().debug() << "Tip " << path.edgesStr() 
				<< " len:" << path.length << " cov:" << path.meanCoverage;

			GraphEdge* targetEdge = path.path.front();
			GraphEdge* complEdge = _graph.complementEdge(targetEdge);

			vecRemove(targetEdge->nodeLeft->outEdges, targetEdge);
			targetEdge->nodeLeft = _graph.addNode();
			targetEdge->nodeLeft->outEdges.push_back(targetEdge);

			//if (targetEdge->selfComplement) continue;

			vecRemove(complEdge->nodeRight->inEdges, complEdge);
			complEdge->nodeRight = _graph.addNode();
			complEdge->nodeRight->inEdges.push_back(complEdge);


			if ((bool)Config::get("remove_alt_edges"))
			{
				for (auto& edge : path.path)
				{
					_graph.removeEdge(_graph.complementEdge(edge));
					_graph.removeEdge(edge);
				}
			}
		}
	}
	outShort = shortClipped;
	outLong = longClipped;
}


