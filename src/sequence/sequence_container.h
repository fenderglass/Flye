//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <string>
#include <limits>

#include "sequence.h"

struct FastaRecord
{
	class Id
	{
	public:
		Id(): _id(std::numeric_limits<uint32_t>::max()) {}

		explicit Id(uint32_t id): _id(id) {}

		bool operator==(const Id& other) const
			{return _id == other._id;}

		bool operator!=(const Id& other) const
			{return !(*this == other);}

		Id rc() const		//reverse complement 
			{return Id(_id + 1 - (_id % 2) * 2);}

		bool strand() const		//true = positive, false = negative
			{return !(_id % 2);}

		size_t hash() const 
		{
			size_t x = _id;
			size_t z = (x += 0x9E3779B97F4A7C15ULL);
			z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
			z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
			return z ^ (z >> 31);
		}

		int signedId() const
			{return (_id % 2) ? -((int)_id + 1) / 2 : (int)_id / 2 + 1;}

		friend std::ostream& operator << (std::ostream& stream, const Id& id)
		{
			stream << std::to_string(id._id);
			return stream;
		}
		
		friend std::istream& operator >> (std::istream& stream, Id& id)
		{
			std::string buffer;
			stream >> buffer;
			id._id = std::stoi(buffer);
			return stream;
		}

		bool operator < (const FastaRecord::Id& other) const
		{
			return _id < other._id;
		}

		friend class SequenceContainer;

	private:
		uint32_t _id;
	};
	static const Id ID_NONE; 
	typedef std::tuple<Id, Id> IdPair;

	FastaRecord(): id(ID_NONE) {}
	FastaRecord(const DnaSequence& sequence, const std::string& description,
				Id id):
		id(id), sequence(sequence), description(description)
	{
	}

	FastaRecord(const FastaRecord& other):
		id(other.id), sequence(other.sequence), 
		description(other.description) {}

	FastaRecord(FastaRecord&& other):
		id (other.id)
	{
		*this = std::move(other);
	}

	FastaRecord& operator=(const FastaRecord& other)
	{
		if (this == &other) return *this;

		id = other.id;
		sequence = other.sequence;
		description = other.description;
		return *this;
	}

	FastaRecord& operator=(FastaRecord&& other)
	{
		if (this == &other) return *this;

		id = other.id;
		sequence = std::move(other.sequence);
		description = std::move(other.description);
		return *this;
	}
	
	Id id;
	DnaSequence sequence;
	std::string description;
};

namespace std
{
	template <>
	struct hash<FastaRecord::Id> 
	{
		size_t operator() (const FastaRecord::Id& h) const throw() 
		{
			 return h.hash();
		}
	};

	template <>
	struct hash<FastaRecord::IdPair> 
	{
		 size_t operator()(const FastaRecord::IdPair& k) const
		 {
			size_t lhs = std::get<0>(k).hash();
			size_t rhs = std::get<1>(k).hash();
			lhs ^= rhs + 0x9ddfea08eb382d69ULL + (lhs << 6) + (lhs >> 2);
			return lhs;
		 }
	};
}

class SequenceContainer
{
public:
	class ParseException : public std::runtime_error 
	{
	public:
		ParseException(const std::string & what):
			std::runtime_error(what)
		{}
	};

	typedef std::vector<FastaRecord> SequenceIndex;

	SequenceContainer():
		_offsetInitialized(false) {}

	void loadFromFile(const std::string& filename, int minReadLength = 0);

	static void writeFasta(const std::vector<FastaRecord>& records,
						   const std::string& fileName,
						   bool  onlyPositiveStrand = false);

	static size_t getMaxSeqId() {return g_nextSeqId;}

	const FastaRecord&  addSequence(const DnaSequence& sequence, 
									const std::string& description);

	const SequenceIndex& iterSeqs() const
	{
		return _seqIndex;
	}

	const FastaRecord& getRecord(FastaRecord::Id seqId) const
	{
		assert(seqId._id - _seqIdOffest < _seqIndex.size());
		assert(_seqIndex[seqId._id - _seqIdOffest].id == seqId);
		return _seqIndex[seqId._id - _seqIdOffest];
	}

	const DnaSequence& getSeq(FastaRecord::Id readId) const
	{
		assert(readId._id - _seqIdOffest < _seqIndex.size());
		assert(_seqIndex[readId._id - _seqIdOffest].id == readId);
		return _seqIndex[readId._id - _seqIdOffest].sequence;
	}

	int32_t seqLen(FastaRecord::Id readId) const
	{
		assert(readId._id - _seqIdOffest < _seqIndex.size());
		assert(_seqIndex[readId._id - _seqIdOffest].id == readId);
		return _seqIndex[readId._id - _seqIdOffest].sequence.length();
	}

	std::string seqName(FastaRecord::Id readId) const
	{
		assert(readId._id - _seqIdOffest < _seqIndex.size());
		assert(_seqIndex[readId._id - _seqIdOffest].id == readId);
		return _seqIndex[readId._id - _seqIdOffest].description;
	}

	int computeNxStat(float fraction) const;

	void   buildPositionIndex();

	size_t globalPosition(FastaRecord::Id seqId, int32_t position) const
	{
		assert(position >= 0 && position < this->seqLen(seqId));
		assert(seqId._id - _seqIdOffest < _seqIndex.size());
		#ifndef NDEBUG
		auto checkGlob = _sequenceOffsets[seqId._id - _seqIdOffest].offset + position;
		FastaRecord::Id checkId;
		int32_t checkPos;
		int32_t outLen;
		this->seqPosition(checkGlob, checkId, checkPos, outLen);
		assert(checkId == seqId && checkPos == position);
		#endif
		return _sequenceOffsets[seqId._id - _seqIdOffest].offset + position;
	}

	const FastaRecord& recordByName(const std::string& name) const
	{
		return this->getRecord(_nameIndex.at(name));
	}

	void seqPosition(size_t globPos, FastaRecord::Id& outSeqId, 
					 int32_t& outPosition, int32_t& outLen) const
	{
		assert(globPos < _sequenceOffsets.back().offset);

		size_t hint = _offsetsHint[globPos / CHUNK];
		while (_sequenceOffsets[hint + 1].offset <= globPos) ++hint;

		outSeqId = FastaRecord::Id(_seqIdOffest + hint);
		outPosition = globPos - _sequenceOffsets[hint].offset;
		outLen = (int32_t)_sequenceOffsets[hint].length;

		assert(outSeqId._id - _seqIdOffest < _seqIndex.size());
		assert(outPosition >= 0 && outPosition < outLen);
		//assert(this->globalPosition(outSeqId, outPosition) == globPos);
	}
	static size_t g_nextSeqId;

private:
	struct OffsetPair
	{
		size_t offset;
		size_t length;
	};

	FastaRecord::Id addSequence(const FastaRecord& sequence);

	size_t readFasta(std::vector<FastaRecord>& record, 
				     const std::string& fileName);

	size_t readFastq(std::vector<FastaRecord>& record, 
				     const std::string& fileName);

	bool   isFasta(const std::string& fileName);

	void   validateSequence(std::string& sequence);

	void   validateHeader(std::string& header);

	SequenceIndex 	_seqIndex;
	size_t 			_seqIdOffest;
	bool   			_offsetInitialized;
	std::unordered_map<std::string, 
					   FastaRecord::Id> _nameIndex;

	//global/local position convertions
	const size_t MAX_SEQUENCE = 1ULL << (8 * 5);
	const size_t CHUNK = 1000;
	std::vector<OffsetPair> _sequenceOffsets;
	std::vector<size_t> 	_offsetsHint;
};

