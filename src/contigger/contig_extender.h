//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "../repeat_graph/repeat_graph.h"
#include "../repeat_graph/read_aligner.h"
#include "../repeat_graph/graph_processing.h"



class ContigExtender
{
public:
	ContigExtender(RepeatGraph& graph, const ReadAligner& aligner,
				   const SequenceContainer& asmSeqs, 
				   const SequenceContainer& readSeqs):
		_graph(graph), _aligner(aligner), 
		_asmSeqs(asmSeqs), _readSeqs(readSeqs) {}

	void generateUnbranchingPaths();
	void generateContigs();
	void outputContigs(const std::string& filename);
	void outputStatsTable(const std::string& filename);
	void outputScaffoldConnections(const std::string& filename);
	//std::vector<UnbranchingPath> getContigPaths();

	const std::vector<UnbranchingPath>& getUnbranchingPaths() 
		{return _unbranchingPaths;}
private:
	struct Contig
	{
		Contig(const UnbranchingPath& corePath):
			graphEdges(corePath), graphPaths({&corePath})
		{}

		UnbranchingPath graphEdges;
		std::vector<const UnbranchingPath*> graphPaths;
		DnaSequence sequence;
	};
	struct Scaffold
	{
		GraphEdge* leftContig;
		GraphEdge* rightContig;
		ska::flat_hash_set<GraphEdge*> repetitiveEdges;
	};
	struct UpathAlignment
	{
		GraphAlignment aln;
		UnbranchingPath* upath;
	};

	std::vector<UnbranchingPath*> asUpaths(const GraphPath& path);
	std::vector<UpathAlignment> asUpathAlignment(const GraphAlignment& aln);

	std::vector<UnbranchingPath> _unbranchingPaths;
	ska::flat_hash_map<GraphEdge*, UnbranchingPath*> _edgeToPath;
	std::vector<Contig> _contigs;

	RepeatGraph& _graph;
	const ReadAligner& _aligner;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
};
