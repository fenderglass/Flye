//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <iostream>
#include <cstring>

#include <cuckoohash_map.hh>

#include "kmer.h"
#include "sequence_container.h"
#include "../common/config.h"
#include "../common/logger.h"


typedef std::map<size_t, size_t> KmerDistribution;

class KmerCounter
{
public:
	KmerCounter(const SequenceContainer& seqContainer):
		_seqContainer(seqContainer), 
		_flatCounter(nullptr), _numKmers(0)
	{}

	~KmerCounter()
	{
		if (_flatCounter) 
		{
			delete[] _flatCounter;
			_flatCounter = nullptr;
		}
	}

	const KmerDistribution& getKmerHist() const
	{
		return _kmerDistribution;
	}

	void   count(bool useFlatCounter);
	size_t getFreq(Kmer kmer) const;
	size_t getKmerNum() const;
	void clear();
	void setOutputProgress(bool progress) {_outputProgress = progress;}

private:
	const SequenceContainer& 	_seqContainer;
	bool _outputProgress;
	bool _useFlatCounter;

	std::atomic<uint8_t>*			_flatCounter;
	//std::vector<std::atomic<char>>  _flatCounter;
	cuckoohash_map<Kmer, size_t> 	_hashCounter;
	KmerDistribution _kmerDistribution;

	std::atomic<size_t> _numKmers;
};

class VertexIndex
{
public:
	~VertexIndex()
	{
		this->clear();
	}
	VertexIndex(const SequenceContainer& seqContainer):
		_seqContainer(seqContainer), _outputProgress(false), 
		_sampleRate(1.0f), _repetitiveFrequency(0),
		_kmerWhitelist(false),
		_kmerCounter(seqContainer)
		//_solidMultiplier(1)
		//_flankRepeatSize(flankRepeatSize)
	{}

	VertexIndex(const VertexIndex&) = delete;
	void operator=(const VertexIndex&) = delete;

private:
	struct IndexChunk
	{
		IndexChunk(): 
			hi(0), low(0) {}
		IndexChunk(const IndexChunk& other): 
			hi(other.hi), low(other.low) {}
		IndexChunk(IndexChunk&& other):
			hi(other.hi), low(other.low) {}
		IndexChunk& operator=(const IndexChunk& other)
		{
			hi = other.hi;
			low = other.low;
			return *this;
		}

		size_t get() const
		{
			return ((size_t)hi << 32) + (size_t)low;
		}
		void set(size_t val)
		{
			low = val & ((1ULL << 32) - 1);
			hi = val >> 32;
		}

		uint8_t hi;
		uint32_t low;
	} __attribute__((packed));
	static_assert(sizeof(IndexChunk) == 5, 
				  "Unexpected size of IndexChunk structure");

	//static const size_t MAX_INDEX = 1ULL << (sizeof(IndexChunk) * 8);

	struct ReadPosition
	{
		ReadPosition(FastaRecord::Id readId = FastaRecord::ID_NONE, 
					 int32_t position = 0):
			readId(readId), position(position) {}
		FastaRecord::Id readId;
		int32_t position;
	};

	struct ReadVector
	{
		ReadVector(uint32_t capacity = 0, uint32_t size = 0): 
			capacity(capacity), size(size), data(nullptr) {}
		uint32_t capacity;
		uint32_t size;
		IndexChunk* data;
	};

public:
	typedef std::map<size_t, size_t> KmerDistribution;

	class KmerPosIterator
	{
	public:
		KmerPosIterator(ReadVector rv, size_t index, bool revComp, 
						const SequenceContainer& seqContainer):
			rv(rv), index(index), revComp(revComp), 
			seqContainer(seqContainer), kmerSize(Parameters::get().kmerSize) 
		{}

		bool operator==(const KmerPosIterator& other) const
		{
			return index == other.index && rv.data == other.rv.data;
		}
		bool operator!=(const KmerPosIterator& other) const
		{
			return !(*this == other);
		}

		//__attribute__((always_inline))
		ReadPosition operator*() const
		{
			size_t globPos = rv.data[index].get();
			FastaRecord::Id seqId;
			int32_t position;
			int32_t seqLen;
			seqContainer.seqPosition(globPos, seqId, position, seqLen);

			if (!revComp)
			{
				return ReadPosition(seqId, position);
			}
			else
			{
				return ReadPosition(seqId.rc(), seqLen - position - kmerSize);
			}
		}

		KmerPosIterator& operator++()
		{
			++index;
			return *this;
		}
	
	private:
		ReadVector rv;
		size_t index;
		bool   revComp;
		const  SequenceContainer& seqContainer;
		size_t kmerSize;
	};

	class IterHelper
	{
	public:
		IterHelper(ReadVector rv, bool revComp, 
				   const SequenceContainer& seqContainer): 
			rv(rv), revComp(revComp), seqContainer(seqContainer) {}

		KmerPosIterator begin()
		{
			return KmerPosIterator(rv, 0, revComp, seqContainer);
		}

		KmerPosIterator end()
		{
			return KmerPosIterator(rv, rv.size, revComp, seqContainer);
		}

	private:
		ReadVector rv;
		bool revComp;
		const SequenceContainer& seqContainer;
	};

	void countKmers();
	void buildIndex(int minCoverage);
	void buildIndexUnevenCoverage(int minCoverage, float selectRate, 
								  int tandemFreq);
	void buildIndexMinimizers(int minCoverage, int wndLen);
	void clear();

	IterHelper iterKmerPos(Kmer kmer) const
	{
		bool revComp = kmer.standardForm();
		return IterHelper(_kmerIndex.find(kmer), revComp,
						  _seqContainer);
	}

	//__attribute__((always_inline))
	/*bool isSolid(Kmer kmer) const
	{
		kmer.standardForm();
		return _kmerIndex.contains(kmer);
	}*/

	bool isRepetitive(Kmer kmer) const
	{
		kmer.standardForm();
		return _repetitiveKmers.contains(kmer);
	}
	
	size_t kmerFreq(Kmer kmer) const
	{
		kmer.standardForm();
		ReadVector rv;
		_kmerIndex.find(kmer, rv);
		return rv.size;
	}

	void outputProgress(bool set) 
	{
		_outputProgress = set;
		_kmerCounter.setOutputProgress(set);
	}

	const KmerDistribution& getKmerHist() const
	{
		return _kmerCounter.getKmerHist();
		//return _kmerDistribution;
	}

	bool isWhitelisted(Kmer kmer) const
	{
		if (!_kmerWhitelist) return true;
		kmer.standardForm();
		return _whitelistedKmers.contains(kmer);
	}

	float getSampleRate() const {return _sampleRate;}

	void loadKmerWhitelist(const std::string& fastaPath);

private:
	//void setRepeatCutoff(int minCoverage);


	struct KmerFreq
	{
		Kmer kmer;
		int32_t position;
		size_t freq;
	};
	std::vector<KmerFreq>
		yieldFrequentKmers(const FastaRecord::Id& seqId,
						   float selctRate, int tandemFreq);

	void allocateIndexMemory();
	void filterFrequentKmers(int minCoverage, float rate);

	const SequenceContainer& _seqContainer;
	//KmerDistribution 		 _kmerDistribution;
	bool    _outputProgress;
	float   _sampleRate;
	size_t  _repetitiveFrequency;
	//int32_t _solidMultiplier;

	const size_t MEM_CHUNK = 32 * 1024 * 1024 / sizeof(IndexChunk);
	std::vector<IndexChunk*> _memoryChunks;

	cuckoohash_map<Kmer, ReadVector> _kmerIndex;
	//cuckoohash_map<Kmer, size_t> 	 _kmerCounts;
	cuckoohash_map<Kmer, char> 	 	 _repetitiveKmers;

	cuckoohash_map<Kmer, char>		 _whitelistedKmers;
	bool 						     _kmerWhitelist;

	KmerCounter _kmerCounter;
};
