//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <limits>
#include <algorithm>
#include <iomanip>
#include <stack>

#include "../common/config.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "extender.h"


Extender::ExtensionInfo Extender::extendContig(FastaRecord::Id startRead)
{

	ska::flat_hash_set<FastaRecord::Id> currentReads;
	currentReads.insert(startRead);
	currentReads.insert(startRead.rc());

	bool rightExtension = true;
	FastaRecord::Id currentRead = startRead;
	std::vector<int> numExtensions;
	std::vector<int> overlapSizes;
	ExtensionInfo exInfo;
	exInfo.reads.push_back(startRead);
	exInfo.assembledLength = _readsContainer.seqLen(startRead);

	auto leftExtendsStart = [startRead, this](const FastaRecord::Id readId)
	{
		for (const auto& ovlp : _ovlpContainer.lazySeqOverlaps(startRead))
		{
			if (ovlp.extId == readId &&
				ovlp.leftShift < -(int)Config::get("maximum_jump")) 
			{
				return true;
			}
		}
		return false;
	};

	while(true)
	{
		const auto& overlaps = _ovlpContainer.lazySeqOverlaps(currentRead);
		std::vector<OverlapRange> extensions;
		//int innerOverlaps = 0;
		for (const auto& ovlp : overlaps)
		{
			//if (_innerReads.contains(ovlp.extId)) ++innerOverlaps;
			if (this->extendsRight(ovlp)) extensions.push_back(ovlp);
		}
		numExtensions.push_back(extensions.size());

		//sort from longes to shortest overlap
		std::sort(extensions.begin(), extensions.end(), 
				  [](const OverlapRange& a, const OverlapRange& b)
					 {return a.curRange() > b.curRange();});

		//checking if read overlaps with one of the used reads
		bool foundExtension = false;

		//getting extension
		int minExtensions = std::max(1, (int)median(numExtensions) / 
										(int)Config::get("max_coverage_drop_rate"));

		/*Logger::get().debug() << _readsContainer.seqName(currentRead) 
			<< "\t" << overlaps.size();
		for (auto& ovlp : overlaps)
		{
			if (this->extendsRight(ovlp))
			{
				Logger::get().debug() << "\t" << 
					_readsContainer.seqName(ovlp.extId) << "\t" 
					<< "o:" << ovlp.curRange() << "\tchim:"
					<< _chimDetector.isChimeric(ovlp.extId)
					<< "\trep:" << this->isRightRepeat(ovlp.extId)
					<< "\text:" << this->countRightExtensions(ovlp.extId);
			}
		}*/

		const OverlapRange* maxExtension = nullptr;
		for (const auto& ovlp : extensions)
		{
			//need to check this condition, otherwise there will
			//be complications when we initiate extension of startRead
			//to the left
			if (leftExtendsStart(ovlp.extId)) continue;
			if (_ovlpContainer.hasSelfOverlaps(ovlp.extId)) continue;

			//try to find a good one
			if (!_chimDetector.isChimeric(ovlp.extId) &&
				this->countRightExtensions(ovlp.extId) >= minExtensions)
			{
				foundExtension = true;
				maxExtension = &ovlp;
				exInfo.assembledLength += ovlp.rightShift;
				currentRead = ovlp.extId;
				break;
			}

			if (!maxExtension || maxExtension->rightShift < ovlp.rightShift)
			{
				maxExtension = &ovlp;
			}
		}
		//in case of suspicious extension make the farthest jump possible
		if (!foundExtension && maxExtension)
		{
			++exInfo.numSuspicious;
			exInfo.assembledLength += maxExtension->rightShift;
			foundExtension = true;
			currentRead = maxExtension->extId;
		}

		//bool overlapsVisited = _innerReads.contains(currentRead);
		if (foundExtension) 
		{
			exInfo.reads.push_back(currentRead);
			overlapSizes.push_back(maxExtension->curRange());
			//overlapsVisited |= currentReads.count(currentRead);
		}
		else
		{
			rightExtension ? exInfo.leftTip = true : exInfo.rightTip = true;
		}

		if (!foundExtension || _innerReads.contains(currentRead) ||
			currentReads.count(currentRead))
		{
			//Logger::get().debug() << "Not found: " << !foundExtension << 
			//	" overlaps visited: " << overlapsVisited;

			//right extension done, try to extend left from start read
			if (rightExtension && !exInfo.reads.empty())
			{
				exInfo.stepsToTurn = exInfo.reads.size();
				rightExtension = false;
				currentRead = exInfo.reads.front().rc();
				std::reverse(exInfo.reads.begin(), exInfo.reads.end());
				for (size_t i = 0; i < exInfo.reads.size(); ++i) 
				{
					exInfo.reads[i] = exInfo.reads[i].rc();
				}
			}
			//done with extension
			else
			{
				break;
			}
		}

		currentReads.insert(currentRead);
		currentReads.insert(currentRead.rc());
	}

	if (!numExtensions.empty())
	{
		exInfo.meanOverlaps = median(numExtensions);
	}
	if (!overlapSizes.empty())
	{
		exInfo.avgOverlapSize = median(overlapSizes);
		exInfo.minOverlapSize = *std::min_element(overlapSizes.begin(), 
												  overlapSizes.end());
	}

	if (_innerReads.contains(exInfo.reads.front())) 
	{
		exInfo.leftAsmOverlap = _readsContainer.seqLen(exInfo.reads.front());
	}
	if (_innerReads.contains(exInfo.reads.back())) 
	{
		exInfo.rightAsmOverlap = _readsContainer.seqLen(exInfo.reads.back());
	}

	return exInfo;
}


void Extender::assembleContigs()
{
	//static const int MAX_JUMP = Config::get("maximum_jump");
	Logger::get().info() << "Extending reads";
	_chimDetector.estimateGlobalCoverage();
	_ovlpContainer.overlapDivergenceStats();
	_innerReads.clear();
	cuckoohash_map<FastaRecord::Id, size_t> coveredReads;
	
	int totalReads = 0;
	for (const auto& read : _readsContainer.iterSeqs())
	{
		if ((int)read.sequence.length() > 
			Parameters::get().minimumOverlap) ++totalReads;
	}
	
	std::mutex indexMutex;
	auto processRead = [this, &indexMutex, &coveredReads, totalReads] 
		(FastaRecord::Id startRead)
	{
		//most of the reads will fall into the inner categoty -
		//so no further processing will be needed
		if (_innerReads.contains(startRead)) return;

		coveredReads.insert(startRead);
		coveredReads.insert(startRead.rc());

		//getting overlaps without caching first - so we don't
		//store overlap information for many trashy reads
		//that won't result into contig extension
		auto startOvlps = _ovlpContainer.quickSeqOverlaps(startRead, 
														  /*max overlaps*/ 100);
		std::vector<OverlapRange> revOverlaps;
		revOverlaps.reserve(startOvlps.size());
		for (const auto& ovlp : startOvlps) revOverlaps.push_back(ovlp.complement());

		int numInnerOvlp = 0;
		int totalOverlaps = 0;
		for (const auto& ovlp : startOvlps)
		{
			if (_innerReads.contains(ovlp.extId)) ++numInnerOvlp;
			++totalOverlaps;
		}

		int maxStartExt = _chimDetector.getOverlapCoverage() * 10;
		//int minStartExt = std::max(1, _chimDetector.getOverlapCoverage() / 50);
		int minStartExt = 1;
		int extLeft = this->countRightExtensions(revOverlaps);
		int extRight = this->countRightExtensions(startOvlps);

		if (_chimDetector.isChimeric(startRead, startOvlps) ||
			std::max(extLeft, extRight) > maxStartExt ||
			std::min(extLeft, extRight) < minStartExt ||
			numInnerOvlp > totalOverlaps / 2) return;
		
		//Good to go!
		ExtensionInfo exInfo = this->extendContig(startRead);

		if (exInfo.reads.size() - exInfo.numSuspicious < 
			(size_t)Config::get("min_reads_in_contig")) return;

		/*if (exInfo.leftAsmOverlap + exInfo.rightAsmOverlap > 
			exInfo.assembledLength + 2 * Parameters::get().minimumOverlap)
		{
			//Logger::get().debug() << "No novel sequence";
			return;
		}*/

		//Exclusive part - updating the overall assembly
		std::lock_guard<std::mutex> guard(indexMutex);
		
		int innerCount = 0;
		//do not count first and last reads - they are inner by defalut
		assert(exInfo.reads.size() >= 4);
		for (size_t i = 1; i < exInfo.reads.size() - 1; ++i)
		{
			if (_innerReads.contains(exInfo.reads[i])) ++innerCount;
		}
		int innerThreshold = std::min((int)Config::get("max_inner_reads"),
									  int((float)Config::get("max_inner_fraction") * 
										  exInfo.reads.size()));
		if (innerCount > innerThreshold)
		{
			Logger::get().debug() << "Discarded contig with "
				<< exInfo.reads.size() << " reads and "
				<< innerCount << " inner overlaps";
			return;
		}

		Logger::get().debug() << "Assembled contig " 
			<< std::to_string(_readLists.size() + 1)
			<< "\n\tWith " << exInfo.reads.size() << " reads"
			<< "\n\tStart read: " << _readsContainer.seqName(startRead)
			<< "\n\tAt position: " << exInfo.stepsToTurn
			<< "\n\tleftTip: " << exInfo.leftTip 
			<< " rightTip: " << exInfo.rightTip
			<< "\n\tSuspicious: " << exInfo.numSuspicious
			<< "\n\tMean extensions: " << exInfo.meanOverlaps
			<< "\n\tAvg overlap len: " << exInfo.avgOverlapSize
			<< "\n\tMin overlap len: " << exInfo.minOverlapSize
			<< "\n\tInner reads: " << innerCount
			<< "\n\tLength: " << exInfo.assembledLength;

		Logger::get().debug() << "Ovlp index size: " << _ovlpContainer.indexSize();
		
		//update inner read index
		ska::flat_hash_set<FastaRecord::Id> rightExtended;
		ska::flat_hash_set<FastaRecord::Id> leftExtended;
		std::vector<OverlapRange> allOverlaps;
		for (const auto& readId : exInfo.reads)
		{
			coveredReads.insert(readId, true);
			coveredReads.insert(readId.rc(), true);
			_innerReads.insert(readId, true);
			_innerReads.insert(readId.rc(), true);

			//repetitive read - will bring to many "off-target" reads
			//int maxExtensions = exInfo.meanOverlaps * 2;
			//if (this->countRightExtensions(readId) > maxExtensions||
			//	this->countRightExtensions(readId.rc()) > maxExtensions) continue;

			//so each read is covered from the left and right
			for (const auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
			{
				allOverlaps.push_back(ovlp);
				coveredReads.insert(ovlp.extId, true);
				coveredReads.insert(ovlp.extId.rc(), true);
			}
		}
		for (const auto& read : this->getInnerReads(allOverlaps))
		{
			_innerReads.insert(read, true);
			_innerReads.insert(read.rc(), true);
		}

		Logger::get().debug() << "Inner: " << 
			_innerReads.size() << " covered: " << coveredReads.size()
			<< " total: "<< totalReads;
		
		_readLists.push_back(std::move(exInfo));
	};

	std::function<void(const FastaRecord::Id&)> threadWorker = 
		[processRead] (const FastaRecord::Id& readId)
	{
		processRead(readId);
	};
	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _readsContainer.iterSeqs())
	{
		if (seq.sequence.length() > (size_t)Parameters::get().minimumOverlap &&
			seq.id.strand())
		{
			allReads.push_back(seq.id);
		}
	}
	std::random_shuffle(allReads.begin(), allReads.end());
	processInParallel(allReads, threadWorker,
					  Parameters::get().numThreads, true);
	//_ovlpContainer.ensureTransitivity(/*only max*/ true);

	bool addSingletons = (bool)Config::get("add_unassembled_reads");
	if (addSingletons)
	{
		std::vector<FastaRecord::Id> sortedByLength;
		for (const auto& seq : _readsContainer.iterSeqs())
		{
			if (seq.id.strand() && !_innerReads.contains(seq.id) &&
				_readsContainer.seqLen(seq.id) > Parameters::get().minimumOverlap)
			{
				sortedByLength.push_back(seq.id);
			}
		}
		std::sort(sortedByLength.begin(), sortedByLength.end(),
				  [this](FastaRecord::Id idOne, FastaRecord::Id idTwo)
				  	{return _readsContainer.seqLen(idOne) > 
							_readsContainer.seqLen(idTwo);});

		int singletonsAdded = 0;
		ska::flat_hash_set<FastaRecord::Id> coveredLocal;
		for (const auto& readId : sortedByLength)
		{
			if (!coveredLocal.count(readId))
			{
				for (const auto& ovlp : _ovlpContainer.lazySeqOverlaps(readId))
				{
					if (ovlp.leftShift >= 0 && ovlp.rightShift <= 0)
					{
						coveredLocal.insert(ovlp.extId);
						coveredLocal.insert(ovlp.extId.rc());
					}
				}
				ExtensionInfo path;
				path.singleton = true;
				path.reads.push_back(readId);
				_readLists.push_back(path);
				++singletonsAdded;
			}
		}
		Logger::get().info() << "Added " << singletonsAdded << " singleton reads";
	}

	this->convertToContigs();
	Logger::get().info() << "Assembled " << _contigPaths.size() << " draft contigs";
}


std::vector<FastaRecord::Id> 
	Extender::getInnerReads(const std::vector<OverlapRange>& ovlps)
{
	static const int WINDOW = Config::get("chimera_window");
	static const int OVERHANG = Config::get("maximum_overhang");

	ska::flat_hash_map<FastaRecord::Id, 
					   std::vector<int32_t>> readsCoverage;
	for (const auto& ovlp: ovlps)
	{
		auto& coverage = readsCoverage[ovlp.extId];
		if (coverage.empty())
		{
			int numWindows = _readsContainer.seqLen(ovlp.extId) / WINDOW;
			if (numWindows < 1)
			{
				Logger::get().warning() << "Wrong read length: " << numWindows;
				numWindows = 1;
			}
			coverage.assign(numWindows, 0);
		}

		for (int pos = ovlp.extBegin / WINDOW + 1; 
			 pos < ovlp.extEnd / WINDOW; ++pos)
		{
			++coverage[pos];
		}
	}

	std::vector<FastaRecord::Id> innerReads;
	for (auto rc : readsCoverage)
	{
		int leftZeros = 0;
		for (size_t i = 0; i < rc.second.size(); ++i)
		{
			if (rc.second[i] != 0) break;
			++leftZeros;
		}
		int rightZeros = 0;
		for (size_t i = 0; i < rc.second.size(); ++i)
		{
			if (rc.second[rc.second.size() - i - 1] != 0) break;
			++rightZeros;
		}
		bool middleZero = false;
		for (size_t i = leftZeros + 1; i < rc.second.size() - rightZeros; ++i)
		{
			if (rc.second[i] == 0) middleZero = true;
		}

		if (!middleZero && leftZeros < OVERHANG / WINDOW && 
			rightZeros < OVERHANG / WINDOW) innerReads.push_back(rc.first);

		/*if (middleZero)
		{
			std::string covStr;
			for (int cov : rc.second)
			{
				covStr += std::to_string(cov) + " ";
			}
			Logger::get().debug() << covStr;
		}*/
	}

	return innerReads;
}

void Extender::convertToContigs()
{
	for (const auto& exInfo : _readLists)
	{
		ContigPath path;
		if (!exInfo.singleton)
		{
			path.name = "contig_" + std::to_string(_contigPaths.size() + 1);
		}
		else
		{
			path.name = "read_" + std::to_string(_contigPaths.size() + 1);
		}
		//path.trimLeft = std::max(0, exInfo.leftAsmOverlap - 
		//							2 * Parameters::get().minimumOverlap);
		//path.trimRight = std::max(0, exInfo.rightAsmOverlap - 
		//							 2 * Parameters::get().minimumOverlap);

		for (size_t i = 0; i < exInfo.reads.size() - 1; ++i)
		{
			bool found = false;
			OverlapRange readsOvlp;

			for (const auto& ovlp : _ovlpContainer.lazySeqOverlaps(exInfo.reads[i]))
			{
				if (ovlp.extId == exInfo.reads[i + 1]) 
				{
					readsOvlp = ovlp;
					found = true;
					break;
				}
			}
			if (!found)
			{
				for (const auto& ovlp : _ovlpContainer.lazySeqOverlaps(exInfo.reads[i + 1]))
				{
					if (ovlp.extId == exInfo.reads[i]) 
					{
						readsOvlp = ovlp.reverse();
						found = true;
						break;
					}
				}
			}
			if (!found) throw std::runtime_error("Ovlp not found!");

			path.sequences.push_back(_readsContainer.getSeq(exInfo.reads[i]));
			path.overlaps.push_back(readsOvlp);
		}
		path.sequences.push_back(_readsContainer.getSeq(exInfo.reads.back()));
		_contigPaths.push_back(std::move(path));
	}
}

int Extender::countRightExtensions(const std::vector<OverlapRange>& ovlps) const
{
	int count = 0;
	for (const auto& ovlp : ovlps)
	{
		if (this->extendsRight(ovlp)) ++count;
	}
	return count;
}

int Extender::countRightExtensions(FastaRecord::Id readId) const
{
	return this->countRightExtensions(_ovlpContainer.lazySeqOverlaps(readId));
}

bool Extender::extendsRight(const OverlapRange& ovlp) const
{
	static const int MAX_JUMP = Config::get("maximum_jump");
	return ovlp.rightShift > MAX_JUMP;
}
