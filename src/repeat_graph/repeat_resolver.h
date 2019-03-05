//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

//Does the "conventional" resolution of bridged repeats.
//Also, classifies edges into unique and repetitive based
//on read alignment

#pragma once

#include "repeat_graph.h"
#include "read_aligner.h"
#include "multiplicity_inferer.h"

class RepeatResolver
{
public:
	RepeatResolver(RepeatGraph& graph, const SequenceContainer& asmSeqs,
				   const SequenceContainer& readSeqs, 
				   ReadAligner& aligner,
				   MultiplicityInferer& multInf): 
		_graph(graph), _asmSeqs(asmSeqs), _readSeqs(readSeqs), 
		_aligner(aligner), _multInf(multInf) {}

	void findRepeats();
	void resolveRepeats();
	void finalizeGraph();

private:
	struct ReadSequence
	{
		FastaRecord::Id readId;
		int32_t start;
		int32_t end;

		int32_t length() const {return end - start;}
	};
	struct Connection
	{
		GraphPath path;
		ReadSequence readSeq;
		int32_t flankLength;
	};

	bool checkByReadExtension(const GraphEdge* edge,
							  const std::vector<GraphAlignment>& alignments);
	void clearResolvedRepeats();
	std::vector<Connection> getConnections();
	int  resolveConnections(const std::vector<Connection>& conns, 
							float minSupport);
	void separatePath(const GraphPath& path, EdgeSequence segment,
					  FastaRecord::Id startId);

	RepeatGraph& _graph;
	const SequenceContainer&   _asmSeqs;
	const SequenceContainer&   _readSeqs;
	ReadAligner& _aligner;
	MultiplicityInferer& _multInf;
	ska::flat_hash_map<GraphEdge*, int> _substractedCoverage;
};
