//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

//Given contigs, represented as list of possibly overlapping sequences
//generate a single consensus sequence for each contig

#pragma once

#include <vector>

#include "../sequence/overlap.h"


struct ContigPath
{
	ContigPath(): trimLeft(0), trimRight(0) {}

	std::string name;
	std::vector<DnaSequence> sequences;
	std::vector<OverlapRange> overlaps;
	int32_t trimLeft;
	int32_t trimRight;
};

class ConsensusGenerator
{
public:
	std::vector<FastaRecord> 
		generateConsensuses(const std::vector<ContigPath>& contigs, 
							bool verbose = true);
	
private:
	struct AlignmentInfo
	{
		std::string alnOne;
		std::string alnTwo;

		int32_t startOne;
		int32_t startTwo;
	};
	typedef ska::flat_hash_map<const OverlapRange*, 
							   AlignmentInfo> AlignmentsMap;

	FastaRecord generateLinear(const ContigPath& path, 
							   const AlignmentsMap& alnMap);
	AlignmentsMap generateAlignments(const std::vector<ContigPath>& contigs, 
									 bool verbose);
	std::pair<int32_t, int32_t> getSwitchPositions(const AlignmentInfo& aln,
												   int32_t prevSwitch);
};
