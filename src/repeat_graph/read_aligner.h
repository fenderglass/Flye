//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)


//Aligns reads to the graph, also updates alignments
//if the graph changes

#pragma once

#include "repeat_graph.h"

struct EdgeAlignment
{
	OverlapRange overlap;
	GraphEdge* edge;
	//EdgeSequence segment;
};
typedef std::vector<EdgeAlignment> GraphAlignment;

class ReadAligner
{
public:
	ReadAligner(RepeatGraph& graph, const SequenceContainer& readSeqs): 
		_graph(graph), _readSeqs(readSeqs) {}

	void alignReads(const std::string& kmerWhitelistPath);
	void updateAlignments();
	const std::vector<GraphAlignment>& getAlignments() const
		{return _readAlignments;}

	void storeAlignments(const std::string& filename);
	void loadAlignments(const std::string& filename);

	typedef std::unordered_map<GraphEdge*, 
					   		   std::vector<GraphAlignment>> AlnIndex;
	AlnIndex makeAlignmentIndex();

	typedef std::unordered_map<GraphEdge*, 
							   std::unordered_map<GraphEdge*, int>> ConnIndex;
	ConnIndex getEdgeConnectivity() const;

private:
	std::vector<GraphAlignment> 
		chainReadAlignments(const std::vector<EdgeAlignment>& ovlps) const;

	float getChainBaseDivergence(const GraphAlignment& aln, bool realign);

	std::vector<GraphAlignment> _readAlignments;

	RepeatGraph& _graph;
	//const SequenceContainer&   _asmSeqs;
	const SequenceContainer&   _readSeqs;
};
