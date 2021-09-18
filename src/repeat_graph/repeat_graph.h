//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <list>
#include <set>

#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../common/config.h"
#include "../common/utils.h"

struct EdgeSequence
{
	EdgeSequence(FastaRecord::Id edgeSeqId = FastaRecord::ID_NONE,
				 int32_t seqLen = 0):
		edgeSeqId(edgeSeqId), seqLen(seqLen), origSeqId(FastaRecord::ID_NONE),
		origSeqLen(0), origSeqStart(0), origSeqEnd(0) {}

		EdgeSequence(FastaRecord::Id origSeq, int32_t origLen,  
					 int32_t origSeqStart, int32_t origSeqEnd):
		edgeSeqId(FastaRecord::ID_NONE), seqLen(origSeqEnd - origSeqStart),
		origSeqId(origSeq), origSeqLen(origLen),
		origSeqStart(origSeqStart), origSeqEnd(origSeqEnd) {}

	EdgeSequence complement() const
	{
		EdgeSequence other(*this);
		other.edgeSeqId = edgeSeqId.rc();
		if (origSeqId != FastaRecord::ID_NONE) 
		{
			other.origSeqId = origSeqId.rc();
			other.origSeqStart = origSeqLen - origSeqEnd - 1;
			other.origSeqEnd = origSeqLen - origSeqStart - 1;
		}
		return other;
	}

	bool operator==(const EdgeSequence& other)
	{
		return edgeSeqId == other.edgeSeqId && 
			   seqLen == other.seqLen &&
			   origSeqId == other.origSeqId &&
			   origSeqLen == other.origSeqLen &&
			   origSeqStart == other.origSeqStart &&
			   origSeqEnd == other.origSeqEnd;
	}

	void dump(std::ostream& os, const SequenceContainer& edgesSeqs)
	{
		std::string origIdString = 
			origSeqId != FastaRecord::ID_NONE ?
			std::to_string(origSeqId.signedId()) : "*";
		os  << edgesSeqs.seqName(edgeSeqId) << " " << seqLen 
			<< " " << origIdString << " " << origSeqLen << " " << origSeqStart 
			<< " " << origSeqEnd;
	}

	void parse(std::istream& is, const SequenceContainer& edgeSeqs)
	{
		std::string edgeSeqName;
		std::string origSeqName;
		is  >> edgeSeqName >> seqLen >> origSeqName >> origSeqLen 
			>> origSeqStart >> origSeqEnd;
		edgeSeqId = edgeSeqs.recordByName(edgeSeqName).id;
		
		if (origSeqName == "*")
		{
			origSeqId = FastaRecord::ID_NONE;
		}
		else
		{
			size_t unsignedId = llabs(atoll(origSeqName.c_str())) * 2 - 2;
			unsignedId += atoll(origSeqName.c_str()) < 0;
			origSeqId = FastaRecord::Id(unsignedId);
		}
	}

	//index in the repeat graph sequence container
	FastaRecord::Id edgeSeqId;
	int32_t seqLen;
	
	//this information is required during repeat graph construction,
	//but not necessary afterwards. It might be used for 
	//some huristics later on, but it it not guaranteed that all
	//edges will have it
	FastaRecord::Id origSeqId;
	int32_t origSeqLen;
	int32_t origSeqStart;
	int32_t origSeqEnd;
};



struct GraphNode;

struct GraphEdge
{
	GraphEdge(GraphNode* nodeLeft, GraphNode* nodeRight, 
			  FastaRecord::Id edgeId = FastaRecord::ID_NONE):
		nodeLeft(nodeLeft), nodeRight(nodeRight), 
		edgeId(edgeId), repetitive(false), 
		selfComplement(false), resolved(false), 
		altHaplotype(false), altGroupId(-1),
		meanCoverage(0), leftLink(nullptr), 
		rightLink(nullptr) {}

	bool isRepetitive() const 
		{return repetitive;}

	bool isLooped() const 
		{return nodeLeft == nodeRight;}

	bool isRightTerminal() const;

	int32_t length() const
	{
		if (seqSegments.empty()) return 0;

		int64_t sumLen = 0;
		for (auto& seqSeg : seqSegments)
		{
			sumLen += seqSeg.seqLen;
		}
		return sumLen / seqSegments.size();
	}

	std::string edgeDescrLong() const
	{
		return "{" + std::to_string(edgeId.signedId()) + "," + std::to_string(length()) + 
			    "," + std::to_string(meanCoverage) + "}";
	}

	std::string edgeDescr() const
	{
		return std::to_string(edgeId.signedId());
	}

	std::unordered_set<GraphEdge*> adjacentEdges();

	/////////////////////////

	GraphNode* nodeLeft;
	GraphNode* nodeRight;

	FastaRecord::Id edgeId;
	std::vector<EdgeSequence> seqSegments;

	bool repetitive;
	bool selfComplement;
	bool resolved;
	bool altHaplotype;
	int  altGroupId;
	//bool unreliable;
	int  meanCoverage;

	GraphEdge* leftLink;
	GraphEdge* rightLink;
};

struct GraphNode
{
	GraphNode(size_t nodeId): nodeId(nodeId) {}

	bool isBifurcation() const
		{return outEdges.size() != 1 || inEdges.size() != 1;}

	std::unordered_set<GraphNode*> neighbors() const
	{
		std::unordered_set<GraphNode*> result;
		for (auto& edge : inEdges) 
		{
			if (edge->nodeLeft != this) result.insert(edge->nodeLeft);
		}
		for (auto& edge : outEdges) 
		{
			if (edge->nodeRight != this) result.insert(edge->nodeRight);
		}

		return result;
	}

	bool isEnd() const
	{
		int inDegree = 0;
		for (auto& edge : inEdges)
		{
			if (!edge->isLooped()) ++inDegree;
		}
		int outDegree = 0;
		for (auto& edge : outEdges)
		{
			if (!edge->isLooped()) ++outDegree;
		}
		return (inDegree == 1 && outDegree == 0) || 
			   (inDegree == 0 && outDegree == 1);
	}

	bool isTelomere() const
	{
		int numIn = 0;
		int numOut = 0;
		for (auto& edge: inEdges)
		{
			if (!edge->isLooped()) ++numIn;
		}
		for (auto& edge: outEdges)
		{
			if (!edge->isLooped()) ++numOut;
		}
		if ((bool)numIn != (bool)numOut)
		{
			return true;
		}
		return false;
	}


	bool isResolved() const
	{
		int inDegree = 0;
		for (auto& edge : inEdges)
		{
			if (!edge->isLooped()) ++inDegree;
		}
		int outDegree = 0;
		for (auto& edge : outEdges)
		{
			if (!edge->isLooped()) ++outDegree;
		}
		return inDegree == 1 && outDegree == 1;
	}

	std::vector<GraphEdge*> inEdges;
	std::vector<GraphEdge*> outEdges;
	size_t nodeId;
};

typedef std::vector<GraphEdge*> GraphPath;


class RepeatGraph
{
public:
	RepeatGraph(const SequenceContainer& asmSeqs, SequenceContainer* graphSeqs):
		 _nextEdgeId(0), _nextNodeId(0), _asmSeqs(asmSeqs), 
		 _edgeSeqsContainer(graphSeqs)
	{}
	~RepeatGraph();

	void build(const std::string& kmerWhitelistPath);
	void updateEdgeSequences();
	void storeGraph(const std::string& filename);
	void loadGraph(const std::string& filename);

	void validateGraph();

	GraphPath  complementPath(const GraphPath& path) const;
	GraphEdge* complementEdge(GraphEdge* edge) const;
	GraphNode* complementNode(GraphNode* node) const;

	//nodes
	GraphNode* addNode()
	{
		GraphNode* node = new GraphNode(_nextNodeId);
		_graphNodes.insert(node);
		++_nextNodeId;
		return node;
	}

	class IterNodes
	{
	public:
		IterNodes(RepeatGraph& graph): _graph(graph) {}

		std::unordered_set<GraphNode*>::iterator begin() 
			{return _graph._graphNodes.begin();}
		std::unordered_set<GraphNode*>::iterator end() 
			{return _graph._graphNodes.end();}
	
	private:
		RepeatGraph& _graph;
	};
	IterNodes iterNodes() {return IterNodes(*this);}
	//
	
	//edges
	GraphEdge* addEdge(GraphEdge&& edge)
	{
		if (this->getEdge(edge.edgeId))
		{
			throw std::runtime_error("Adding edge with duplicated id");
		}

		GraphEdge* newEdge = new GraphEdge(edge);
		newEdge->nodeLeft->outEdges.push_back(newEdge);
		newEdge->nodeRight->inEdges.push_back(newEdge);
		_sortedEdges.insert(newEdge);
		_idToEdge[newEdge->edgeId] = newEdge;
		if (newEdge->selfComplement)
		{
			_idToEdge[newEdge->edgeId.rc()] = newEdge;
		}
		return newEdge;
	}
	/*bool hasEdge(GraphEdge* edge)
	{
		return _idToEdge.count(edge->edgeId);
		//return _sortedEdges.count(edge);
	}*/
	GraphEdge* getEdge(FastaRecord::Id edgeId)
	{
		if (_idToEdge.count(edgeId)) return _idToEdge[edgeId];
		return nullptr;
	}
	class IterEdges
	{
	public:
		IterEdges(RepeatGraph& graph): _graph(graph) {}

		std::set<GraphEdge*>::iterator begin() 
			{return _graph._sortedEdges.begin();}
		std::set<GraphEdge*>::iterator end() 
			{return _graph._sortedEdges.end();}
	
	private:
		RepeatGraph& _graph;
	};

	IterEdges iterEdges() {return IterEdges(*this);}

	void removeEdge(GraphEdge* edge)
	{
		vecRemove(edge->nodeRight->inEdges, edge);
		vecRemove(edge->nodeLeft->outEdges, edge);
		_sortedEdges.erase(edge);
		_idToEdge.erase(edge->edgeId);
		_deletedEdges.insert(edge);
		//delete edge;
	}

	void removeNode(GraphNode* node)
	{
		std::unordered_set<GraphEdge*> toRemove;
		for (auto& edge : node->outEdges) 
		{
			vecRemove(edge->nodeRight->inEdges, edge);
			toRemove.insert(edge);
		}
		for (auto& edge : node->inEdges) 
		{
			vecRemove(edge->nodeLeft->outEdges, edge);
			toRemove.insert(edge);
		}
		for (auto& edge : toRemove)
		{
			_sortedEdges.erase(edge);
			_idToEdge.erase(edge->edgeId);
			_deletedEdges.insert(edge);
			//delete edge;
		}
		_graphNodes.erase(node);
		_deletedNodes.insert(node);
		//delete node;
	}

	//
	FastaRecord::Id newEdgeId()
	{
		size_t curId = _nextEdgeId;
		_nextEdgeId += 2;
		return FastaRecord::Id(curId);
	}

	const SequenceContainer& edgeSequences() {return *_edgeSeqsContainer;}

	EdgeSequence addEdgeSequence(const DnaSequence& sequence, 
							 	 int32_t start, int32_t length,
							 	 const std::string& description);

	void disconnectRight(GraphEdge* edge)
	{
		GraphNode* newNode = this->addNode();
		vecRemove(edge->nodeRight->inEdges, edge);
		edge->nodeRight = newNode;
		edge->nodeRight->inEdges.push_back(edge);
	};

	void disconnectLeft(GraphEdge* edge)
	{
		GraphNode* newNode = this->addNode();
		vecRemove(edge->nodeLeft->outEdges, edge);
		edge->nodeLeft = newNode;
		edge->nodeLeft->outEdges.push_back(edge);
	};

	void linkEdges(GraphEdge* leftEdge, GraphEdge* rightEdge)
	{
		if (leftEdge->rightLink || rightEdge->leftLink)
		{
			if (leftEdge->rightLink != rightEdge ||
				rightEdge->leftLink != leftEdge) Logger::get().warning() << "Relinking!";

			return;
		}

		leftEdge->rightLink = rightEdge;
		rightEdge->leftLink = leftEdge;
	}

private:
	size_t _nextEdgeId;
	size_t _nextNodeId;

	struct GluePoint
	{
		GluePoint(size_t id = 0, FastaRecord::Id seqId = FastaRecord::ID_NONE,
				  int32_t position = 0):
			pointId(id), seqId(seqId), position(position) {}

		size_t 	pointId;
		FastaRecord::Id seqId;
		int32_t	position;
	};

	struct RepeatCluster
	{
		RepeatCluster(FastaRecord::Id seqId = FastaRecord::ID_NONE,
					  size_t clusterId = 0, int32_t start = 0,
					  int32_t end = 0):
			seqId(seqId), clusterId(clusterId), start(start), end(end) {}

		FastaRecord::Id seqId;
		size_t clusterId;
		int32_t start;
		int32_t end;
	};

	void getGluepoints(OverlapContainer& ovlps);
	void initializeEdges(const OverlapContainer& asmOverlaps);
	void collapseTandems();
	void logEdges();
	void checkGluepointProjections(const OverlapContainer& asmOverlaps);
	
	const SequenceContainer& _asmSeqs;
	SequenceContainer* 		 _edgeSeqsContainer;
	const int _maxSeparation = Config::get("max_separation");

	std::unordered_map<FastaRecord::Id, 
					   std::vector<GluePoint>> _gluePoints;

	std::unordered_map<FastaRecord::Id, GraphEdge*> _idToEdge;
	std::unordered_set<GraphEdge*> _deletedEdges;
	std::unordered_set<GraphNode*> _graphNodes;
	std::unordered_set<GraphNode*> _deletedNodes;

	struct CmpId 
	{
		bool operator() (GraphEdge* const &e1, GraphEdge* const &e2) const
		{return e1->edgeId < e2->edgeId;}
	};
	std::set<GraphEdge*, CmpId> _sortedEdges;
};
