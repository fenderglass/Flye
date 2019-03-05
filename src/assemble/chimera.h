//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "../sequence/overlap.h"
#include "../sequence/sequence_container.h"

class ChimeraDetector
{
public:
	ChimeraDetector(const SequenceContainer& readContainer,
					OverlapContainer& ovlpContainer):
		_seqContainer(readContainer),
		_ovlpContainer(ovlpContainer), 
		_overlapCoverage(0)
	{}

	void estimateGlobalCoverage();
	bool isChimeric(FastaRecord::Id readId);
	bool isChimeric(FastaRecord::Id readId, 
					const std::vector<OverlapRange>& readOvlps);
	int  getOverlapCoverage() const {return _overlapCoverage;}
	int  getRightTrim(FastaRecord::Id readId);

private:
	std::vector<int32_t> 
		getReadCoverage(FastaRecord::Id readId,
						const std::vector<OverlapRange>& readOvlps);
	bool testReadByCoverage(FastaRecord::Id readId,
							const std::vector<OverlapRange>& readOvlps);

	const SequenceContainer& _seqContainer;
	OverlapContainer& _ovlpContainer;

	cuckoohash_map<FastaRecord::Id, bool> _chimeras;
	int _overlapCoverage;
};
