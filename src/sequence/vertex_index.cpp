//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <cmath>

#include "vertex_index.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/config.h"
#include "../common/memory_info.h"


void VertexIndex::countKmers()
{
	_kmerCounter.count(/*use flat counter*/ true);
}


void VertexIndex::buildIndexUnevenCoverage(int globalMinFreq, float selectRate,
										   int tandemFreq)
{
	//this->setRepeatCutoff(globalMinFreq);

	//_solidMultiplier = 1;

	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}

	//first, count the number of k-mers that will be actually stored in the index
	_kmerIndex.reserve(_kmerCounter.getKmerNum() / 10);
	if (_outputProgress) Logger::get().info() << "Filling index table (1/2)";
	std::function<void(const FastaRecord::Id&)> initializeIndex = 
	[this, globalMinFreq, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto topKmers = this->yieldFrequentKmers(readId, selectRate, tandemFreq);
		for (auto kmerFreq : topKmers)
		{
			kmerFreq.kmer.standardForm();

			if (kmerFreq.freq < (size_t)globalMinFreq) continue;

			ReadVector defVec((uint32_t)1, (uint32_t)0);
			_kmerIndex.upsert(kmerFreq.kmer, 
							  [](ReadVector& rv){++rv.capacity;}, defVec);
		}
	};
	processInParallel(allReads, initializeIndex, 
					  Parameters::get().numThreads, _outputProgress);
	
	this->filterFrequentKmers(globalMinFreq, (float)Config::get("repeat_kmer_rate"));
	this->allocateIndexMemory();

	if (_outputProgress) Logger::get().info() << "Filling index table (2/2)";
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, globalMinFreq, selectRate, tandemFreq] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto topKmers = this->yieldFrequentKmers(readId, selectRate, tandemFreq);
		for (auto kmerFreq : topKmers)
		{
			if (kmerFreq.freq < (size_t)globalMinFreq ||
				kmerFreq.freq > _repetitiveFrequency) continue;

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

			//will not trigger update for k-mer not in the index
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

	_kmerCounter.clear();

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
	Logger::get().debug() << "Mean k-mer index frequency: " 
		<< (float)totalEntries / _kmerIndex.size();
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

void VertexIndex::filterFrequentKmers(int minCoverage, float rate)
{
	size_t totalKmers = 0;
	size_t uniqueKmers = 0;
	for (const auto& kmer : _kmerIndex.lock_table())
	{
		if (kmer.second.capacity >= (size_t)minCoverage)
		{
			totalKmers += kmer.second.capacity;
			uniqueKmers += 1;
		}
	}
	float meanFrequency = (float)totalKmers / (uniqueKmers + 1);
	_repetitiveFrequency = rate * meanFrequency;
	
	size_t repetitiveKmers = 0;
	for (const auto& kmer : _kmerIndex.lock_table())
	{
		if (kmer.second.capacity > _repetitiveFrequency)
		{
			//++repetitiveKmers;
			repetitiveKmers += kmer.second.capacity;
			_repetitiveKmers.insert(kmer.first, true);
		}
	}

	for (const auto& kmer : _repetitiveKmers.lock_table())
	{
		_kmerIndex.erase(kmer.first);
	}

	float filteredRate = (float)repetitiveKmers / totalKmers;
	Logger::get().debug() << "Mean k-mer frequency: " 
						  << meanFrequency;
	Logger::get().debug() << "Repetitive k-mer frequency: " 
						  << _repetitiveFrequency;
	Logger::get().debug() << "Filtered " << repetitiveKmers 
						  << " repetitive k-mers (" <<
						  filteredRate << ")";
}

std::vector<VertexIndex::KmerFreq>
	VertexIndex::yieldFrequentKmers(const FastaRecord::Id& seqId,
									float selectRate, int tandemFreq)
{
	thread_local std::unordered_map<Kmer, size_t> localFreq;
	localFreq.clear();
	std::vector<KmerFreq> topKmers;
	topKmers.reserve(_seqContainer.seqLen(seqId));

	for (const auto& kmerPos : IterKmers(_seqContainer.getSeq(seqId)))
	{
		auto stdKmer = kmerPos.kmer;
		stdKmer.standardForm();
		size_t freq = _kmerCounter.getFreq(stdKmer);

		++localFreq[stdKmer];
		topKmers.push_back({kmerPos.kmer, kmerPos.position, freq});
	}

	if (topKmers.empty()) return {};
	std::sort(topKmers.begin(), topKmers.end(),
			  [](const KmerFreq& k1, const KmerFreq& k2)
			   {return k1.freq > k2.freq;});
	const size_t maxKmers = selectRate * topKmers.size();
	const size_t minFreq = topKmers[maxKmers].freq;

	auto itVec = topKmers.begin();
	while(itVec != topKmers.end() && itVec->freq >= minFreq) ++itVec;
	topKmers.erase(itVec, topKmers.end());

	if (tandemFreq > 0)
	{
		topKmers.erase(std::remove_if(topKmers.begin(), topKmers.end(),
							[tandemFreq](KmerFreq kf)
							{
								kf.kmer.standardForm();
								return localFreq[kf.kmer] > (size_t)tandemFreq;
							}), 
					   topKmers.end());
	}

	return topKmers;
}



void VertexIndex::allocateIndexMemory()
{
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
}

void VertexIndex::buildIndexMinimizers(int minCoverage, int wndLen)
{
	if (_outputProgress) Logger::get().info() << "Building minimizer index";

	std::vector<FastaRecord::Id> allReads;
	size_t totalLen = 0;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
		if (seq.id.strand()) totalLen += seq.sequence.length();
	}

	_kmerIndex.reserve(1000000);
	if (_outputProgress) Logger::get().info() << "Pre-calculating index storage";
	std::function<void(const FastaRecord::Id&)> initializeIndex = 
	[this, minCoverage, wndLen] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;

		auto minimizers = yieldMinimizers(_seqContainer.getSeq(readId), wndLen);
		for (auto kmerPos : minimizers)
		{
			auto stdKmer = kmerPos.kmer;
			stdKmer.standardForm();
			ReadVector defVec((uint32_t)1, (uint32_t)0);
			_kmerIndex.upsert(stdKmer, 
							  [](ReadVector& rv){++rv.capacity;}, defVec);
		}
	};
	processInParallel(allReads, initializeIndex, 
					  Parameters::get().numThreads, _outputProgress);

	this->filterFrequentKmers(minCoverage, (float)Config::get("repeat_kmer_rate"));
	this->allocateIndexMemory();
	
	if (_outputProgress) Logger::get().info() << "Filling index";
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, minCoverage, wndLen] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;
		auto minimizers = yieldMinimizers(_seqContainer.getSeq(readId), wndLen);
		for (auto kmerPos : minimizers)
		{
			FastaRecord::Id targetRead = readId;
			bool revCmp = kmerPos.kmer.standardForm();
			if (revCmp)
			{
				kmerPos.position = _seqContainer.seqLen(readId) - 
										kmerPos.position -
										Parameters::get().kmerSize;
				targetRead = targetRead.rc();
			}

			if (_repetitiveKmers.contains(kmerPos.kmer)) continue;

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
	Logger::get().debug() << "K-mer index size: " << totalEntries;
	Logger::get().debug() << "Mean k-mer frequency: " 
		<< (float)totalEntries / _kmerIndex.size();

	float minimizerRate = (float)totalLen / totalEntries;
	Logger::get().debug() << "Minimizer rate: " << minimizerRate;
	_sampleRate = minimizerRate;
}


void VertexIndex::clear()
{
	for (auto& chunk : _memoryChunks) delete[] chunk;
	_memoryChunks.clear();

	_kmerIndex.clear();
	_kmerIndex.reserve(0);

	_kmerCounter.clear();
	//_kmerCounts.reserve(0);
}


void KmerCounter::count(bool useFlatCounter)
{
	//Logger::get().debug() << "Before counter: " 
	//	<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	if (useFlatCounter && Parameters::get().kmerSize > 17)
	{
		throw std::runtime_error("Can't use flat counter for k-mer size > 17");
	}
	_useFlatCounter = useFlatCounter;

	//flat array for all possible k-mers, 4 bits for each
	//in case of k=17, takes 8Gb
	static const size_t COUNTER_LEN = std::pow(4, Parameters::get().kmerSize) / 2;
	if (useFlatCounter)
	{
		_flatCounter = new std::atomic<uint8_t>[COUNTER_LEN];
		std::memset(_flatCounter, 0, COUNTER_LEN);
	}
 
	if (_outputProgress) Logger::get().info() << "Counting k-mers:";
	std::function<void(const FastaRecord::Id&)> readUpdate = 
	[this] (const FastaRecord::Id& readId)
	{
		if (!readId.strand()) return;
		
		for (auto kmerPos : IterKmers(_seqContainer.getSeq(readId)))
		{
			kmerPos.kmer.standardForm();
			bool addOne = true;
			if (_useFlatCounter)
			{
				size_t arrayPos = kmerPos.kmer.numRepr() / 2;
				bool highBits = kmerPos.kmer.numRepr() % 2;

				while (true)
				{
					uint8_t expected = _flatCounter[arrayPos]; 
					uint8_t count = highBits ? (expected >> 4) : (expected & 15);
					if (count == 15)
					{
						break;
					}

					uint8_t updated = highBits ? (expected + 16) : (expected + 1);
					if (_flatCounter[arrayPos].compare_exchange_weak(expected,  updated))
					{
						if (count == 0) ++_numKmers;
						addOne = false; //not saturated yet, don't update hash counter
						break;
					}
				}
			}

			if (addOne)
			{
				_hashCounter.upsert(kmerPos.kmer, [](size_t& num){++num;}, 1);
			}
		}
	};
	std::vector<FastaRecord::Id> allReads;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		allReads.push_back(seq.id);
	}
	processInParallel(allReads, readUpdate, Parameters::get().numThreads, _outputProgress);

	Logger::get().debug() << "Updating k-mer histogram";
	if (_useFlatCounter)
	{
		for (size_t kmerId = 0; kmerId < COUNTER_LEN * 2; ++kmerId)
		{
			Kmer kmer(kmerId);
			size_t freq = this->getFreq(kmer);
			if (freq > 0) _kmerDistribution[freq] += 1;
			//if (kmerId % 1000000 == 0) Logger::get().debug() << kmerId << " " << freq;
		}
	}
	else
	{
		for (const auto& kmer : _hashCounter.lock_table())
		{
			_kmerDistribution[kmer.second] += 1;
		}
	}

	//Logger::get().debug() << "After counter: " 
	//	<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	Logger::get().debug() << "Hash size: " << _hashCounter.size();
	Logger::get().debug() << "Total k-mers " << _numKmers;
}


size_t KmerCounter::getFreq(Kmer kmer) const
{
	//kmer.standardForm();

	size_t addCount = 0;
	if (_useFlatCounter)
	{
		size_t arrayPos = kmer.numRepr() / 2;
		bool highBits = kmer.numRepr() % 2;
		uint8_t count = highBits ? (_flatCounter[arrayPos]) >> 4 : (_flatCounter[arrayPos] & 15);
		if (count < 15) 
		{
			return count;
		}
		else
		{
			addCount = count;
		}
	}

	size_t freq = 0;
	_hashCounter.find(kmer, freq);
	return freq + addCount;
}

void KmerCounter::clear()
{
	_hashCounter.clear();
	_hashCounter.reserve(0);
	if (_flatCounter)
	{
		delete[] _flatCounter;
		_flatCounter = nullptr;
	}
}

size_t KmerCounter::getKmerNum() const
{
	if (!_useFlatCounter) return _hashCounter.size();
	return _numKmers;
}


void VertexIndex::loadKmerWhitelist(const std::string& fastaPath)
{
	std::vector<FastaRecord> kmers;
	SequenceContainer s;
	s.readFasta(kmers, fastaPath);
	for (auto& rec : kmers)
	{
		if (rec.sequence.length() != Parameters::get().kmerSize)
		{
			throw std::runtime_error("kmer size mismatch");
		}
		Kmer kmer(rec.sequence, 0, Parameters::get().kmerSize);
		kmer.standardForm();
		_whitelistedKmers.insert(kmer, true);
	}
	_kmerWhitelist = true;
}
