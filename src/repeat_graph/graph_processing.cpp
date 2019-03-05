//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <deque>

#include "graph_processing.h"
#include "../common/logger.h"
#include "../common/config.h"
#include "../common/utils.h"

//a helper function that calls
//all other simplification procedures
void GraphProcessor::simplify()
{
	//this->trimTips();
	this->fixChimericJunctions();
	this->condenceEdges();
	this->collapseBulges();
	this->condenceEdges();
	//this->trimTips();
}

//finds and removes graph structures that
//originate from the typical chimeric read pattern 
//(when a read contain two consecivtive reversed copies of the real sequence)
void GraphProcessor::fixChimericJunctions()
{
	//a very specific case: 1 in - 1 out
	ska::flat_hash_set<GraphNode*> simpleCases;
	for (auto& node : _graph.iterNodes())
	{
		if (!node->isBifurcation() &&
			(node->inEdges.front()->edgeId ==
			 node->outEdges.front()->edgeId.rc()))
		{
			simpleCases.insert(node);
		}
	}
	for (auto& node : simpleCases)
	{
		GraphNode* newNode = _graph.addNode();
		GraphEdge* cutEdge = node->outEdges.front();
		newNode->outEdges.push_back(cutEdge);
		cutEdge->nodeLeft = newNode;
		node->outEdges.clear();
	}

	//more common case: 2 in - 2 out
	ska::flat_hash_set<GraphNode*> complexCases;
	for (auto& node : _graph.iterNodes())
	{
		if (node->inEdges.size() != 2 ||
			node->outEdges.size() != 2) continue;
		auto& inEdges = node->inEdges;
		auto& outEdges = node->outEdges;
		if (inEdges[0]->edgeId.rc() != outEdges[0]->edgeId)
		{
			//match INs with OUTs
			std::swap(inEdges[0], inEdges[1]);
		}

		if (inEdges[0]->edgeId.rc() == outEdges[0]->edgeId &&
			inEdges[1]->edgeId.rc() == outEdges[1]->edgeId)
		{
			complexCases.insert(node);
		}
	}
	for(auto& node : complexCases)
	{
		GraphNode* newNode = _graph.addNode();
		node->inEdges[1]->nodeRight = newNode;
		node->outEdges[0]->nodeLeft = newNode;
		newNode->inEdges.push_back(node->inEdges[1]);
		newNode->outEdges.push_back(node->outEdges[0]);

		node->inEdges.pop_back();
		node->outEdges.erase(node->outEdges.begin());
	}

	Logger::get().debug() << "Removed " 
		<< simpleCases.size() << " simple and " << complexCases.size()
		<< " double chimeric junctions";
}

//Collapses simple small bulges
void GraphProcessor::collapseBulges()
{
	const int MAX_BUBBLE = Parameters::get().minimumOverlap;
	ska::flat_hash_set<std::pair<GraphNode*, GraphNode*>,
					   pairhash> toFix;
	//finding the bubbles
	for (auto& edge : _graph.iterEdges())
	{
		if (edge->isLooped()) continue;
		std::vector<GraphEdge*> parallelEdges;
		for (auto& parEdge : edge->nodeLeft->outEdges)
		{
			if (parEdge->nodeRight == edge->nodeRight) 
			{
				parallelEdges.push_back(parEdge);
			}
		}
		if (parallelEdges.size() != 2) continue;
		if (parallelEdges[0]->edgeId == parallelEdges[1]->edgeId.rc()) continue;
		if (parallelEdges[0]->length() > MAX_BUBBLE || 
			parallelEdges[1]->length() > MAX_BUBBLE) continue;

		toFix.emplace(edge->nodeLeft, edge->nodeRight);
	}

	//collapsing them by leaving only one branch and
	//adding the other's sequences as alternative to the
	//remaining branch
	for (auto& nodes : toFix)
	{
		std::vector<GraphEdge*> parallelEdges;
		for (auto& parEdge : nodes.first->outEdges)
		{
			if (parEdge->nodeRight == nodes.second) 
			{
				parallelEdges.push_back(parEdge);
			}
		}
		GraphEdge* edgeOne = parallelEdges[0];
		GraphEdge* edgeTwo = parallelEdges[1];
		if (abs(edgeOne->edgeId.signedId()) > abs(edgeTwo->edgeId.signedId()))
		{
			std::swap(edgeOne, edgeTwo);
		}
		for (auto& seg : edgeTwo->seqSegments)
		{
			edgeOne->seqSegments.push_back(seg);
		}
		_graph.removeEdge(edgeTwo);
	}
	Logger::get().debug() << "Collapsed " << toFix.size() / 2 << " bulges";
}

//Removing tips
/*void GraphProcessor::trimTips()
{
	const int TIP_THRESHOLD = Config::get("tip_length_threshold");
	ska::flat_hash_set<GraphEdge*> toRemove;
	for (GraphEdge* tipEdge : _graph.iterEdges())
	{
		int rightOut = tipEdge->nodeRight->outEdges.size();
		int leftIn = tipEdge->nodeLeft->inEdges.size();
		int leftOut = tipEdge->nodeLeft->outEdges.size();

		if (tipEdge->length() < TIP_THRESHOLD && rightOut == 0 &&
			leftIn == 1 && leftOut == 2 &&
			!tipEdge->nodeLeft->inEdges.front()->isLooped())
		{
			toRemove.insert(tipEdge);
		}
	}

	for (auto& edge : toRemove)
	{
		vecRemove(edge->nodeLeft->outEdges, edge);
		edge->nodeLeft = _graph.addNode();
		edge->nodeLeft->outEdges.push_back(edge);

		GraphEdge* complEdge = _graph.complementEdge(edge);
		vecRemove(complEdge->nodeRight->inEdges, complEdge);
		complEdge->nodeRight = _graph.addNode();
		complEdge->nodeRight->inEdges.push_back(complEdge);
	}
	Logger::get().debug() << toRemove.size() << " tips clipped";
}*/

//This function collapses non-branching edges paths in the graph.
//The tricky part is the sequence representation of the new edges.
//Two (or more) consecutive egdes will only be collapsed into one
//if there exist at least one consigous sub-sequence from the input assembly
//that can represent this edge
void GraphProcessor::condenceEdges()
{
	int edgesRemoved = 0;
	int edgesAdded = 0;

	//helper function that checks if the simplification is
	//possible and returns the new edges that should replace
	//the original ones
	auto collapseEdges = [] (const GraphPath& edges)
	{
		std::vector<GraphEdge> newEdges;
		std::list<EdgeSequence> growingSeqs(edges.front()->seqSegments.begin(),
											edges.front()->seqSegments.end());
		assert(edges.size() > 1);
		size_t prevStart = 0;
		for (size_t i = 1; i < edges.size(); ++i)
		{
			auto prevSeqs = growingSeqs;
			for (auto prevSeg = growingSeqs.begin(); 
				 prevSeg != growingSeqs.end(); )
			{
				bool continued = false;
				for (auto& nextSeg : edges[i]->seqSegments)
				{
					if (prevSeg->origSeqId == nextSeg.origSeqId &&
						prevSeg->origSeqEnd == nextSeg.origSeqStart)
					{
						continued = true;
						prevSeg->origSeqEnd = nextSeg.origSeqEnd;
					}
				}
				if (!continued)
				{
					prevSeg = growingSeqs.erase(prevSeg);
				}
				else
				{
					++prevSeg;
				}
			}

			if (growingSeqs.empty())
			{
				newEdges.emplace_back(edges[prevStart]->nodeLeft, 
									  edges[i - 1]->nodeRight);
				std::copy(prevSeqs.begin(), prevSeqs.end(),
				  		  std::back_inserter(newEdges.back().seqSegments));

				std::copy(edges[i]->seqSegments.begin(), 
						  edges[i]->seqSegments.end(), 
						  std::back_inserter(growingSeqs));
				prevStart = i;
			}
		}

		newEdges.emplace_back(edges[prevStart]->nodeLeft, 
							  edges.back()->nodeRight);
		std::copy(growingSeqs.begin(), growingSeqs.end(),
				  std::back_inserter(newEdges.back().seqSegments));

		return newEdges;
	};

	//Try to collapse each unbranching path. 
	//If it is possible - replace with new edges
	auto toCollapse = this->getUnbranchingPaths();
	for (auto& unbranchingPath : toCollapse)
	{
		if (!unbranchingPath.id.strand()) continue;
		if (unbranchingPath.path.size() == 1) continue;

		GraphPath complPath = _graph.complementPath(unbranchingPath.path);
		auto newEdges = collapseEdges(unbranchingPath.path);
		if (newEdges.size() == unbranchingPath.path.size()) continue;

		std::string addedStr;
		for (auto& edge : newEdges)
		{
			//do not add collapsed short loops
			//if (edge.length() < Parameters::get().minimumOverlap &&
			//	edge.isLooped()) continue;	

			GraphEdge addFwd = edge;
			addFwd.edgeId = _graph.newEdgeId();

			GraphEdge addRev(_graph.complementNode(edge.nodeRight),
							 _graph.complementNode(edge.nodeLeft),
							 addFwd.edgeId.rc());

			//complementary segments
			for (auto seqSeg : addFwd.seqSegments)
			{
				addRev.seqSegments.push_back(seqSeg.complement());
			}

			GraphEdge* addedEdge = _graph.addEdge(std::move(addFwd));
			_graph.addEdge(std::move(addRev));

			addedStr += std::to_string(addedEdge->edgeId.signedId()) + " -> ";
		}

		if (addedStr.size() > 4) addedStr.erase(addedStr.size() - 4);
		//Logger::get().debug() << "Collapsed: " << unbranchingPath.edgesStr() 
		//	<< " to " << addedStr;

		ska::flat_hash_set<GraphEdge*> toRemove;
		for (auto& edge : unbranchingPath.path) toRemove.insert(edge);
		for (auto& edge : complPath) toRemove.insert(edge);
		for (auto& edge : toRemove) _graph.removeEdge(edge);

		edgesRemoved += unbranchingPath.path.size();
		edgesAdded += newEdges.size();
	}

	//Logger::get().debug() << "Removed " << edgesRemoved << " edges";
	//Logger::get().debug() << "Added " << edgesAdded << " edges";
	Logger::get().debug() << "Collapsed " << edgesRemoved - edgesAdded << " edges";
}

//converts edges to unbranching paths
std::vector<UnbranchingPath> GraphProcessor::getEdgesPaths() const
{
	std::vector<UnbranchingPath> paths;
	for (auto& edge : _graph.iterEdges())
	{
		GraphPath path = {edge};
		paths.emplace_back(path, edge->edgeId, false,
						   edge->length(), edge->meanCoverage);
		paths.back().repetitive = edge->repetitive;
	}
	return paths;
}

//Finds unbranching paths
std::vector<UnbranchingPath> GraphProcessor::getUnbranchingPaths() const
{
	ska::flat_hash_map<FastaRecord::Id, size_t> edgeIds;
	size_t nextEdgeId = 0;
	auto pathToId = [&edgeIds, &nextEdgeId](GraphPath path)
	{
		if (!edgeIds.count(path.front()->edgeId))
		{
			for (auto edge : path)
			{
				edgeIds[edge->edgeId] = nextEdgeId;
				edgeIds[edge->edgeId.rc()] = nextEdgeId + 1;
			}
			nextEdgeId += 2;
			return FastaRecord::Id(nextEdgeId - 2);
		}
		return FastaRecord::Id(edgeIds[path.front()->edgeId]);
	};
	
	std::vector<UnbranchingPath> unbranchingPaths;
	ska::flat_hash_set<GraphEdge*> visitedEdges;
	for (auto edge : _graph.iterEdges())
	{
		if (visitedEdges.count(edge)) continue;
		visitedEdges.insert(edge);

		GraphPath traversed;
		traversed.push_back(edge);
		if (!edge->selfComplement)
		{
			GraphNode* curNode = edge->nodeLeft;
			while (!curNode->isBifurcation() &&
				   !curNode->inEdges.empty() &&
				   !visitedEdges.count(curNode->inEdges.front()) &&
				   !curNode->inEdges.front()->selfComplement)
			{
				traversed.push_back(curNode->inEdges.front());
				visitedEdges.insert(traversed.back());
				curNode = curNode->inEdges.front()->nodeLeft;
			}
			std::reverse(traversed.begin(), traversed.end());
			curNode = edge->nodeRight;
			while (!curNode->isBifurcation() &&
				   !curNode->outEdges.empty() &&
				   !visitedEdges.count(curNode->outEdges.front()) &&
				   !curNode->outEdges.front()->selfComplement)
			{
				traversed.push_back(curNode->outEdges.front());
				visitedEdges.insert(traversed.back());
				curNode = curNode->outEdges.front()->nodeRight;
			}
		}

		FastaRecord::Id edgeId = pathToId(traversed);
		bool circular = (traversed.front()->nodeLeft == 
							traversed.back()->nodeRight) &&
						traversed.front()->nodeLeft->outEdges.size() == 1 &&
						traversed.front()->nodeLeft->inEdges.size() == 1;

		bool repetitive = traversed.front()->isRepetitive() || 
						  traversed.back()->isRepetitive();

		int64_t contigLength = 0;
		int64_t sumCov = 0;
		for (auto& edge : traversed) 
		{
			contigLength += edge->length();
			sumCov += (int64_t)edge->meanCoverage * (int64_t)edge->length();
		}
		int meanCoverage = contigLength ? sumCov / contigLength : 0;

		unbranchingPaths.emplace_back(traversed, edgeId, circular, 
							  contigLength, meanCoverage);
		unbranchingPaths.back().repetitive = repetitive;
	}
	return unbranchingPaths;
}
