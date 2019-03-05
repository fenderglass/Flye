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
	//const int SHORT_EDGE = Config::get("unique_edge_length");

	//alternative coverage
	ska::flat_hash_map<GraphEdge*, std::vector<int32_t>> wndCoverage;

	for (auto& edge : _graph.iterEdges())
	{
		size_t numWindows = edge->length() / WINDOW;
		wndCoverage[edge].assign(numWindows, 0);
	}

	for (auto& path : _aligner.getAlignments())
	{
		for (size_t i = 0; i < path.size(); ++i)
		{
			auto& ovlp = path[i].overlap;
			auto& coverage = wndCoverage[path[i].edge];
			for (int pos = ovlp.extBegin / WINDOW + 1; 
			 	 pos < ovlp.extEnd / WINDOW; ++pos)
			{
				if (pos >= 0 && 
					pos < (int)coverage.size())
				{
					++coverage[pos];
				}
			}
		}
	}

	int64_t sumCov = 0;
	int64_t sumLength = 0;
	for (auto& edgeCoverage : wndCoverage)
	{
		//if (edgeCoverage.first->length() < SHORT_EDGE) continue;
		for (auto& cov : edgeCoverage.second)
		{
			sumCov += (int64_t)cov;
			++sumLength;
		}
	}
	_meanCoverage = (sumLength != 0) ? sumCov / sumLength : 1;

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

	_uniqueCovThreshold = 2;
	if (!edgesCoverage.empty())
	{
		const float MULT = 1.75f;	//at least 1.75x of mean coverage
		_uniqueCovThreshold = MULT * quantile(edgesCoverage, 75);
	}
	Logger::get().debug() << "Unique coverage threshold " << _uniqueCovThreshold;
}

//removes edges with low coverage support from the graph
void MultiplicityInferer::removeUnsupportedEdges()
{
	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	int32_t coverageThreshold = std::round((float)this->getMeanCoverage() / 
											Config::get("graph_cov_drop_rate"));
	coverageThreshold = std::max(1, coverageThreshold);
	if (Parameters::get().unevenCoverage)
	{
		coverageThreshold = std::min(coverageThreshold, 2);
	}
	Logger::get().debug() << "Read coverage cutoff: " << coverageThreshold;

	ska::flat_hash_set<GraphEdge*> edgesRemove;
	for (auto& path : unbranchingPaths)
	{
		if (!path.id.strand()) continue;

		if (path.meanCoverage <= coverageThreshold)
		{
			Logger::get().debug() << "Low coverage: " 
				<< path.edgesStr() << " " << path.meanCoverage;
			for (auto& edge : path.path)
			{
				edgesRemove.insert(edge);
				edgesRemove.insert(_graph.complementEdge(edge));
			}
		}
	}
	for (auto& edge : edgesRemove) _graph.removeEdge(edge);
	Logger::get().debug() << "Removed " << edgesRemove.size() / 2
		<< " unsupported edges";

	_aligner.updateAlignments();
}

void MultiplicityInferer::removeUnsupportedConnections()
{
	ska::flat_hash_map<GraphEdge*, int32_t> rightConnections;
	ska::flat_hash_map<GraphEdge*, int32_t> leftConnections;

	for (auto& readPath : _aligner.getAlignments())
	{
		if (readPath.size() < 2) continue;
		//int overhang = std::max(readPath.front().overlap.curBegin,
		//						readPath.back().overlap.curLen - 
		//							readPath.back().overlap.curEnd);
		//if (overhang > (int)Config::get("maximum_overhang")) continue;

		for (size_t i = 0; i < readPath.size() - 1; ++i)
		{
			if (readPath[i].edge == readPath[i + 1].edge &&
				readPath[i].edge->isLooped()) continue;
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

	auto disconnectRight = [this](GraphEdge* edge)
	{
		GraphNode* newNode = _graph.addNode();
		vecRemove(edge->nodeRight->inEdges, edge);
		edge->nodeRight = newNode;
		edge->nodeRight->inEdges.push_back(edge);
	};
	auto disconnectLeft = [this](GraphEdge* edge)
	{
		GraphNode* newNode = _graph.addNode();
		vecRemove(edge->nodeLeft->outEdges, edge);
		edge->nodeLeft = newNode;
		edge->nodeLeft->outEdges.push_back(edge);
	};

	for (auto& edge : _graph.iterEdges())
	{
		if (!edge->edgeId.strand() || edge->isLooped()) continue;
		GraphEdge* complEdge = _graph.complementEdge(edge);

		int32_t coverageThreshold = edge->meanCoverage / 
								Config::get("graph_cov_drop_rate");

		//Logger::get().debug() << "Adjacencies: " << edge->edgeId.signedId() << " "
		//	<< leftConnections[edge] / 2 << " " << rightConnections[edge] / 2;

		if (!edge->nodeRight->isEnd() &&
			edge->nodeRight->isBifurcation() &&
			rightConnections[edge] / 2 < coverageThreshold)
		{
			Logger::get().debug() << "Chimeric right: " <<
				edge->edgeId.signedId() << " " << rightConnections[edge] / 2;

			disconnectRight(edge);
			disconnectLeft(complEdge);

			if (edge->selfComplement) continue;	//already discinnected
		}
		if (!edge->nodeLeft->isEnd() &&
			edge->nodeLeft->isBifurcation() &&
			leftConnections[edge] / 2 < coverageThreshold)
		{
			Logger::get().debug() << "Chimeric left: " <<
				edge->edgeId.signedId() << " " << leftConnections[edge] / 2;

			disconnectLeft(edge);
			disconnectRight(complEdge);
		}
	}

	//GraphProcessor proc(_graph, _asmSeqs);
	//proc.trimTips();
	_aligner.updateAlignments();
}

//collapse loops (that also could be viewed as
//bubbles with one branch of length 0)
void MultiplicityInferer::collapseHeterozygousLoops()
{
	const float COV_MULT = 1.5;

	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	ska::flat_hash_set<FastaRecord::Id> toUnroll;
	ska::flat_hash_set<FastaRecord::Id> toRemove;
	for (auto& loop : unbranchingPaths)
	{
		if (!loop.isLooped()) continue;
		if (loop.path.front()->selfComplement) continue;

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
		if (loop.meanCoverage > COV_MULT * entrancePath->meanCoverage ||
			loop.meanCoverage > COV_MULT * entrancePath->meanCoverage) continue;

		//loop should not be longer than other branches
		if (loop.length > entrancePath->length ||
			loop.length > exitPath->length) continue;

		if (loop.meanCoverage < 
			(entrancePath->meanCoverage + exitPath->meanCoverage) / 4)
		{
			toRemove.insert(loop.id);
			toRemove.insert(loop.id.rc());
		}
		else
		{
			toUnroll.insert(loop.id);
			toUnroll.insert(loop.id.rc());
		}
	}

	for (auto& path : unbranchingPaths)
	{
		if (toUnroll.count(path.id))
		{
			GraphNode* newNode = _graph.addNode();
			size_t id = path.nodeLeft()->inEdges[0] == path.path.front();
			GraphEdge* prevEdge = path.nodeLeft()->inEdges[id];

			vecRemove(path.nodeLeft()->outEdges, path.path.front());
			vecRemove(path.nodeLeft()->inEdges, prevEdge);
			path.nodeLeft() = newNode;
			newNode->outEdges.push_back(path.path.front());
			prevEdge->nodeRight = newNode;
			newNode->inEdges.push_back(prevEdge);
		}
		if (toRemove.count(path.id))
		{
			GraphNode* newLeft = _graph.addNode();
			GraphNode* newRight = _graph.addNode();

			vecRemove(path.nodeLeft()->outEdges, path.path.front());
			vecRemove(path.nodeLeft()->inEdges, path.path.back());
			path.nodeLeft() = newLeft;
			newRight->inEdges.push_back(path.path.back());
			path.nodeRight() = newRight;
			newLeft->outEdges.push_back(path.path.front());
		}
	}

	Logger::get().debug() << "Unrolled " << toUnroll.size() / 2
		<< " heterozygous loops";
	Logger::get().debug() << "Removed " << toRemove.size() / 2
		<< " heterozygous loops";
	_aligner.updateAlignments();
}

void MultiplicityInferer::trimTips()
{
	const int MAX_TIP = 100000;
	const int MIN_COV_DIFF = 3;

	//const int TIP_THRESHOLD = Config::get("tip_length_threshold");
	ska::flat_hash_set<FastaRecord::Id> toRemove;
	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();
	ska::flat_hash_map<GraphEdge*, UnbranchingPath*> ubIndex;
	for (auto& path : unbranchingPaths)
	{
		for (auto& edge: path.path) ubIndex[edge] = &path;
	}

	for (auto& tipPath : unbranchingPaths)
	{
		if (tipPath.nodeRight()->outEdges.size() > 0) continue;
		if (tipPath.length > MAX_TIP) continue;

		GraphNode* tipNode = tipPath.nodeLeft();

		//get "non-tip" node's entrance and exit
		std::vector<UnbranchingPath*> entrances;
		for (GraphEdge* edge : tipNode->inEdges)
		{
			//for (auto& path : unbranchingPaths)
			//{
			UnbranchingPath& path = *ubIndex[edge];
			if (path.path.back() == edge && !path.isLooped() &&
				path.nodeLeft()->inEdges.size() > 0) entrances.push_back(&path);
			//}
		}
		std::vector<UnbranchingPath*> exits;
		for (GraphEdge* edge : tipNode->outEdges)
		{
			//for (auto& path : unbranchingPaths)
			//{
			UnbranchingPath& path = *ubIndex[edge];
			if (path.path.front() == edge && !path.isLooped() &&
				path.nodeRight()->outEdges.size() > 0) exits.push_back(&path);
			//}
		}
		if (entrances.size() != 1 || exits.size() != 1) continue;

		if (entrances.front()->meanCoverage > MIN_COV_DIFF * tipPath.meanCoverage &&
			exits.front()->meanCoverage > MIN_COV_DIFF * tipPath.meanCoverage)
		{
			toRemove.insert(tipPath.id);
		}
	}

	for (auto& path : unbranchingPaths)
	{
		if (toRemove.count(path.id))
		{
			GraphEdge* targetEdge = path.path.front();
			GraphEdge* complEdge = _graph.complementEdge(targetEdge);

			vecRemove(targetEdge->nodeLeft->outEdges, targetEdge);
			targetEdge->nodeLeft = _graph.addNode();
			targetEdge->nodeLeft->outEdges.push_back(targetEdge);

			vecRemove(complEdge->nodeRight->inEdges, complEdge);
			complEdge->nodeRight = _graph.addNode();
			complEdge->nodeRight->inEdges.push_back(complEdge);
		}
	}
	Logger::get().debug() << toRemove.size() << " tips clipped";
	_aligner.updateAlignments();
}


void MultiplicityInferer::collapseHeterozygousBulges()
{
	const float MAX_COV_VAR = 0.20;
	const float MAX_LEN_VAR = 0.50;

	GraphProcessor proc(_graph, _asmSeqs);
	auto unbranchingPaths = proc.getUnbranchingPaths();

	ska::flat_hash_set<FastaRecord::Id> toSeparate;
	for (auto& path : unbranchingPaths)
	{
		if (path.isLooped()) continue;

		std::vector<UnbranchingPath*> twoPaths;
		for (auto& candEdge : unbranchingPaths)
		{
			if (candEdge.nodeLeft() == path.nodeLeft() &&
				candEdge.nodeRight() == path.nodeRight()) 
			{
				twoPaths.push_back(&candEdge);
			}
		}

		//making sure the structure is ok
		if (twoPaths.size() != 2) continue;
		if (twoPaths[0]->id == twoPaths[1]->id.rc()) continue;
		if (toSeparate.count(twoPaths[0]->id) || 
			toSeparate.count(twoPaths[1]->id)) continue;
		if (twoPaths[0]->nodeLeft()->inEdges.size() != 1 ||
			twoPaths[0]->nodeRight()->outEdges.size() != 1) continue;

		UnbranchingPath* entrancePath = nullptr;
		UnbranchingPath* exitPath = nullptr;
		for (auto& cand : unbranchingPaths)
		{
			if (cand.nodeRight() == 
				twoPaths[0]->nodeLeft()) entrancePath = &cand;
			if (cand.nodeLeft() == twoPaths[0]->nodeRight()) exitPath = &cand;
		}

		//coverage requirement: sum over two branches roughly equals to
		//exit and entrance coverage
		float covSum = twoPaths[0]->meanCoverage + twoPaths[1]->meanCoverage;
		float entranceDiff = fabsf(covSum - entrancePath->meanCoverage) / covSum;
		float exitDiff = fabsf(covSum - exitPath->meanCoverage) / covSum;
		if (entranceDiff > MAX_COV_VAR || exitDiff > MAX_COV_VAR) continue;

		//length requirement: branches have roughly the same length
		//and are significantly shorter than entrance/exits
		if (abs(twoPaths[0]->length - twoPaths[1]->length) >
			MAX_LEN_VAR * std::min(twoPaths[0]->length, 
					 			   twoPaths[1]->length)) continue;
		float bubbleSize = (twoPaths[0]->length + twoPaths[1]->length) / 2;
		if (bubbleSize > entrancePath->length ||
			bubbleSize > exitPath->length) continue;

		if (twoPaths[0]->meanCoverage < twoPaths[1]->meanCoverage)
		{
			toSeparate.insert(twoPaths[0]->id);
			toSeparate.insert(twoPaths[0]->id.rc());
		}
		else
		{
			toSeparate.insert(twoPaths[1]->id);
			toSeparate.insert(twoPaths[1]->id.rc());
		}
	}

	for (auto& path : unbranchingPaths)
	{
		if (toSeparate.count(path.id))
		{
			GraphNode* newLeft = _graph.addNode();
			GraphNode* newRight = _graph.addNode();
			vecRemove(path.nodeLeft()->outEdges, path.path.front());
			vecRemove(path.nodeRight()->inEdges, path.path.back());
			path.nodeLeft() = newLeft;
			path.nodeRight() = newRight;
			newLeft->outEdges.push_back(path.path.front());
			newRight->inEdges.push_back(path.path.back());
		}
	}

	Logger::get().debug() << "Popped " << toSeparate.size() / 2 
		<< " heterozygous bulges";
	_aligner.updateAlignments();
}
