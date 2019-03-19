//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <queue>

#include "vertex_index.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/config.h"


void VertexIndex::countKmers(size_t hardThreshold, int genomeSize)
{
	if (Parameters::get().kmerSize > 31)
	{
		throw std::runtime_error("Maximum k-mer size is 31");
	}

	Logger::get().debug() << "Hard threshold set to " << hardThreshold;
	if (hardThreshold == 0)
	{
		throw std::runtime_error("Wrong hard threshold value: " + 
								 std::to_string(hardThreshold));
	}
	Logger::get().debug() << "Started k-mer counting";

	size_t preCountSize = 1024 * 1024 * 1024;	//1G by default
	if (genomeSize > (int)Config::get("big_genome_threshold"))
	{
		preCountSize *= 4 * 4;					//16G in case of larger genomes
	}
	auto preCounters = new std::atomic<unsigned char>[preCountSize];
	for (size_t i = 0; i < preCountSize; ++i) preCounters[i] = 0;

	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}

	//first pass: filling up naive hash counting filter
	if (_outputProgress) Logger::get().info() << "Counting k-mers (1/2):";
	std::function<void(const FastaRecord::Id&)> preCountUpdate = 
	[&preCounters, hardThreshold, this, preCountSize] 
		(const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;
		
		int32_t nextKmerPos = _sampleRate;
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			if (_sampleRate > 1) //subsampling
			{
				if (--nextKmerPos > 0) continue;
				nextKmerPos = _sampleRate + 
					(int32_t)((kmerPos.kmer.hash() ^ readId.hash()) % 3) - 1;
			}

			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
			}
			size_t kmerBucket = kmerPos.kmer.hash() % preCountSize;

			unsigned char expected = 0;
			while (true)
			{
				expected = preCounters[kmerBucket]; 
				if (expected == std::numeric_limits<unsigned char>::max()) 
				{
					break;
				}
				if (preCounters[kmerBucket]
						.compare_exchange_weak(expected, expected + 1))
				{
					break;
				}
			}
		}
	};
	processInParallel(allReads, preCountUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	//second pass: counting kmers that have passed the filter
	if (_outputProgress) Logger::get().info() << "Counting k-mers (2/2):";

	std::function<void(const FastaRecord::Id&)> countUpdate = 
	[&preCounters, hardThreshold, this, preCountSize] 
		(const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		int32_t nextKmerPos = _sampleRate;
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			if (_sampleRate > 1) //subsampling
			{
				if (--nextKmerPos > 0) continue;
				nextKmerPos = _sampleRate + 
					(int32_t)((kmerPos.kmer.hash() ^ readId.hash()) % 3) - 1;
			}

			kmerPos.kmer.standardForm();
			/*if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
			}*/

			size_t count = preCounters[kmerPos.kmer.hash() % preCountSize];
			if (count >= hardThreshold)
			{
				_kmerCounts.upsert(kmerPos.kmer, [](size_t& num){++num;}, 1);
			}
		}
	};
	processInParallel(allReads, countUpdate, 
					  Parameters::get().numThreads, _outputProgress);
	
	for (const auto& kmer : _kmerCounts.lock_table())
	{
		_kmerDistribution[kmer.second] += 1;
		_repetitiveFrequency = std::max(_repetitiveFrequency, kmer.second);
	}
	delete[] preCounters;
}

namespace
{
	struct KmerFreq
	{
		Kmer kmer;
		size_t position;
		size_t freq;
	};
}

void VertexIndex::buildIndexUnevenCoverage(int minCoverage, float selectRate,
										   int tandemFreq)
{
	//_solidMultiplier = 2;
	_solidMultiplier = 1;

	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}

	//first, count the number of k-mers that will be actually stored in the index
	_kmerIndex.reserve(_kmerCounts.size() / 10);
	if (_outputProgress) Logger::get().info() << "Filling index table (1/2)";
	std::function<void(const FastaRecord::Id&)> initializeIndex = 
	[this, minCoverage, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		thread_local std::unordered_map<Kmer, size_t> localFreq;
		localFreq.clear();
		std::vector<KmerFreq> topKmers;
		topKmers.reserve(_seqContainer.seqLen(readId));

		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			kmerPos.kmer.standardForm();
			size_t freq = 1;
			_kmerCounts.find(kmerPos.kmer, freq);

			topKmers.push_back({kmerPos.kmer, (size_t)kmerPos.position, freq});
			++localFreq[kmerPos.kmer];
		}

		if (topKmers.empty()) return;
		std::sort(topKmers.begin(), topKmers.end(),
				  [](const KmerFreq& k1, const KmerFreq& k2)
				   {return k1.freq > k2.freq;});
		const size_t maxKmers = selectRate * topKmers.size();
		const size_t minFreq = std::max((size_t)minCoverage, topKmers[maxKmers].freq);

		for (auto kmerFreq : topKmers)
		{
			if (kmerFreq.freq < minFreq) break;
			if (kmerFreq.freq > _repetitiveFrequency ||
				localFreq[kmerFreq.kmer] > (size_t)tandemFreq) continue;

			ReadVector defVec((uint32_t)1, (uint32_t)0);
			_kmerIndex.upsert(kmerFreq.kmer, 
							  [](ReadVector& rv){++rv.capacity;}, defVec);
		}
	};
	processInParallel(allReads, initializeIndex, 
					  Parameters::get().numThreads, _outputProgress);
	
	_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
	size_t chunkOffset = 0;
	//Important: since packed structures are apparently not thread-safe,
	//make sure that adjacent k-mer index arrays (that are accessed in parallel)
	//do not overlap within 8-byte window
	const size_t PADDING = 1;
	for (auto& kmer : _kmerIndex.lock_table())
	{
		if (MEM_CHUNK < kmer.second.capacity + PADDING) 
		{
			throw std::runtime_error("k-mer is too frequent");
		}
		if (MEM_CHUNK - chunkOffset < kmer.second.capacity + PADDING)
		{
			_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
			chunkOffset = 0;
		}
		kmer.second.data = _memoryChunks.back() + chunkOffset;
		chunkOffset += kmer.second.capacity + PADDING;
	}

	if (_outputProgress) Logger::get().info() << "Filling index table (2/2)";
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, minCoverage, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		thread_local std::unordered_map<Kmer, size_t> localFreq;
		localFreq.clear();
		std::vector<KmerFreq> topKmers;
		topKmers.reserve(_seqContainer.seqLen(readId));

		for (const auto& kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			auto stdKmer = kmerPos.kmer;
			stdKmer.standardForm();
			size_t freq = 1;
			_kmerCounts.find(stdKmer, freq);

			++localFreq[stdKmer];
			topKmers.push_back({kmerPos.kmer, (size_t)kmerPos.position, freq});
		}

		if (topKmers.empty()) return;
		std::sort(topKmers.begin(), topKmers.end(),
				  [](const KmerFreq& k1, const KmerFreq& k2)
				   {return k1.freq > k2.freq;});
		const size_t maxKmers = selectRate * topKmers.size();
		const size_t minFreq = std::max((size_t)minCoverage, topKmers[maxKmers].freq);

		for (auto kmerFreq : topKmers)
		{
			if (kmerFreq.freq < minFreq) break;

			KmerPosition kmerPos(kmerFreq.kmer, kmerFreq.position);
			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}

			if (kmerFreq.freq > _repetitiveFrequency ||
				localFreq[kmerPos.kmer] > (size_t)tandemFreq) continue;

			//in case kmer not in index yet, creates a new vector
			//with a single element in it
			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, this](ReadVector& rv)
				{
					if (rv.size == rv.capacity) 
					{
						Logger::get().warning() << "Index size mismatch " << rv.capacity;
						return;
					}
					size_t globPos = _seqContainer
							.globalPosition(targetRead, kmerPos.position);
					rv.data[rv.size].set(globPos);
					++rv.size;
				});
		}
	};
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	Logger::get().debug() << "Sorting k-mer index";
	for (const auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const IndexChunk& p1, const IndexChunk& p2)
				  	{return p1.get() < p2.get();});
	}
	
	size_t totalEntries = 0;
	for (const auto& kmerRec : _kmerIndex.lock_table())
	{
		totalEntries += kmerRec.second.size;
	}
	Logger::get().debug() << "Selected k-mers: " << _kmerIndex.size();
	Logger::get().debug() << "Index size: " << totalEntries;
}

namespace
{
	template <class T>
	size_t getFreq(T& histIter)
		{return histIter->first;};

	template <class T>
	size_t getCount(T& histIter)
		{return histIter->second;};

}

void VertexIndex::setRepeatCutoff(int minCoverage)
{
	size_t totalKmers = 0;
	size_t uniqueKmers = 0;
	for (auto mapPair = this->getKmerHist().begin();
		 mapPair != this->getKmerHist().end(); ++mapPair)
	{
		if (minCoverage <= (int)getFreq(mapPair))
		{
			totalKmers += getCount(mapPair) * getFreq(mapPair);
			uniqueKmers += getCount(mapPair);
		}
	}
	float meanFrequency = (float)totalKmers / (uniqueKmers + 1);
	_repetitiveFrequency = (float)Config::get("repeat_kmer_rate") * meanFrequency;
	
	size_t repetitiveKmers = 0;
	for (auto mapPair = this->getKmerHist().rbegin();
		 mapPair != this->getKmerHist().rend(); ++mapPair)
	{
		if (getFreq(mapPair) > _repetitiveFrequency)
		{
			repetitiveKmers += getCount(mapPair);
		}
	}
	float filteredRate = (float)repetitiveKmers / uniqueKmers;
	Logger::get().debug() << "Repetitive k-mer frequency: " 
						  << _repetitiveFrequency;
	Logger::get().debug() << "Filtered " << repetitiveKmers 
						  << " repetitive k-mers (" <<
						  filteredRate << ")";
}

void VertexIndex::buildIndex(int minCoverage)
{
	if (_outputProgress) Logger::get().info() << "Filling index table";
	_solidMultiplier = 1;
	
	//"Replacing" k-mer couns with k-mer index. We need multiple passes
	//to avoid peaks in memory usage during the hash table extensions +
	//prevent memory fragmentation
	
	size_t kmerEntries = 0;
	size_t solidKmers = 0;
	for (const auto& kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second &&
			kmer.second < _repetitiveFrequency)
		{
			kmerEntries += kmer.second;
			++solidKmers;
		}
		if (kmer.second > _repetitiveFrequency)
		{
			_repetitiveKmers.insert(kmer.first, true);
		}
	}
	Logger::get().debug() << "Sampling rate: " << _sampleRate;
	Logger::get().debug() << "Solid k-mers: " << solidKmers;
	Logger::get().debug() << "K-mer index size: " << kmerEntries;
	Logger::get().debug() << "Mean k-mer frequency: " 
		<< (float)kmerEntries / solidKmers;

	_kmerIndex.reserve(solidKmers);
	for (const auto& kmer : _kmerCounts.lock_table())
	{
		if ((size_t)minCoverage <= kmer.second &&
			kmer.second < _repetitiveFrequency)
		{
			ReadVector rv((uint32_t)kmer.second, 0);
			_kmerIndex.insert(kmer.first, rv);
		}
	}
	_kmerCounts.clear();
	_kmerCounts.reserve(0);

	_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
	size_t chunkOffset = 0;
	//Important: since packed structures are apparently not thread-safe,
	//make sure that adjacent k-mer index arrays (that are accessed in parallel)
	//do not overlap within 8-byte window
	const size_t PADDING = 1;
	for (auto& kmer : _kmerIndex.lock_table())
	{
		if (MEM_CHUNK < kmer.second.capacity + PADDING) 
		{
			throw std::runtime_error("k-mer is too frequent");
		}
		if (MEM_CHUNK - chunkOffset < kmer.second.capacity + PADDING)
		{
			_memoryChunks.push_back(new IndexChunk[MEM_CHUNK]);
			chunkOffset = 0;
		}
		kmer.second.data = _memoryChunks.back() + chunkOffset;
		chunkOffset += kmer.second.capacity + PADDING;
	}
	//Logger::get().debug() << "Total chunks " << _memoryChunks.size()
	//	<< " wasted space: " << wasted;

	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		int32_t nextKmerPos = _sampleRate;
		//int32_t seqLen = _seqContainer.seqLen(readId);
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			if (_sampleRate > 1) //subsampling
			{
				if (--nextKmerPos > 0) continue;
				nextKmerPos = _sampleRate + 
					(int32_t)((kmerPos.kmer.hash() ^ readId.hash()) % 3) - 1;
			}

			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}
			
			_kmerIndex.update_fn(kmerPos.kmer, 
				[targetRead, &kmerPos, this](ReadVector& rv)
				{
					size_t globPos = _seqContainer
							.globalPosition(targetRead, kmerPos.position);
					//if (globPos > MAX_INDEX) throw std::runtime_error("Too much!");
					rv.data[rv.size].set(globPos);
					//rv.data[rv.size] = ReadPosition(targetRead, 
					//								kmerPos.position);
					++rv.size;
				});
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}
	processInParallel(allReads, indexUpdate, 
					  Parameters::get().numThreads, _outputProgress);

	Logger::get().debug() << "Sorting k-mer index";
	for (const auto& kmerVec : _kmerIndex.lock_table())
	{
		std::sort(kmerVec.second.data, kmerVec.second.data + kmerVec.second.size,
				  [](const IndexChunk& p1, const IndexChunk& p2)
				  	{return p1.get() < p2.get();});
	}
}


void VertexIndex::clear()
{
	for (auto& chunk : _memoryChunks) delete[] chunk;
	_memoryChunks.clear();

	_kmerIndex.clear();
	_kmerIndex.reserve(0);

	_kmerCounts.clear();
	_kmerCounts.reserve(0);
}
