//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <memory>

#include "sequence_container.h"
#include "../common/config.h"

static_assert(sizeof(size_t) == 8, "32-bit architectures are not supported");

class Kmer
{
public:
	Kmer(): _representation(0) {}

	Kmer(const DnaSequence& dnaString, 
		   size_t start, size_t length):
		_representation(0)
	{
		if (length != Parameters::get().kmerSize)
		{
			throw std::runtime_error("Kmer length inconsistency");
		}

		for (size_t i = start; i < start + length; ++i)	
		{
			_representation <<= 2;
			_representation += dnaString.atRaw(i);
		}
	}

	Kmer reverseComplement()
	{
		KmerRepr tmpRepr = _representation;
		Kmer newKmer;

		for (unsigned int i = 0; i < Parameters::get().kmerSize; ++i)
		{
			newKmer._representation <<= 2;
			newKmer._representation += ~tmpRepr & 3;
			tmpRepr >>= 2;
		}

		return newKmer;
	}

	bool standardForm()
	{
		Kmer complKmer = this->reverseComplement();
		if (complKmer._representation < _representation)
		{
			_representation = complKmer._representation;
			return true;
		}
		return false;
	}

	void appendRight(DnaSequence::NuclType dnaSymbol)
	{
		_representation <<= 2;
		_representation += dnaSymbol;

		KmerRepr kmerSize = Parameters::get().kmerSize;
		KmerRepr kmerMask = ((KmerRepr)1 << kmerSize * 2) - 1;
		_representation &= kmerMask;
	}

	void appendLeft(DnaSequence::NuclType dnaSymbol)
	{
		_representation >>= 2;

		KmerRepr kmerSize = Parameters::get().kmerSize;
		KmerRepr shift = kmerSize * 2 - 2;
		_representation += dnaSymbol << shift;
	}

	typedef size_t KmerRepr;

	bool operator == (const Kmer& other) const
		{return this->_representation == other._representation;}
	bool operator != (const Kmer& other) const
		{return !(*this == other);}
	size_t hash() const
	{
		size_t x = _representation;
		size_t z = (x += 0x9E3779B97F4A7C15ULL);
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
		return z ^ (z >> 31);
	}
	bool operator< (const Kmer& other)
	{
		return _representation < other._representation;
	}

private:
	KmerRepr _representation;
};

namespace std
{
	template <>
	struct hash<Kmer>
	{
		std::size_t operator()(const Kmer& kmer) const
		{
			return kmer.hash();
		}
	};
}

struct KmerPosition
{
	KmerPosition(Kmer kmer, int32_t position):
		kmer(kmer), position(position) {}
	Kmer kmer;
	int32_t position;
};

class KmerIterator
{
public:
    typedef std::forward_iterator_tag iterator_category;

	KmerIterator(const DnaSequence* readSeq, size_t position):
		_readSeq(readSeq),
		_position(position)
	{
		if (position != readSeq->length() - Parameters::get().kmerSize)
		{
			//_kmer = Kmer(readSeq->substr(0, Parameters::get().kmerSize));
			_kmer = Kmer(*readSeq, 0, Parameters::get().kmerSize);
		}
	}

	bool operator==(const KmerIterator& other) const
	{
		return _readSeq == other._readSeq && _position == other._position;
	}

	bool operator!=(const KmerIterator& other) const
	{
		return !(*this == other);
	}

	KmerPosition operator*() const
	{
		return KmerPosition(_kmer, _position);
	}

	KmerIterator& operator++()
	{
		size_t appendPos = _position + Parameters::get().kmerSize;
		_kmer.appendRight(_readSeq->atRaw(appendPos));
		++_position;
		return *this;
	}

protected:
	const DnaSequence* _readSeq;
	size_t 	_position;
	Kmer 	_kmer;
};


class IterKmers
{
public:
	IterKmers(const DnaSequence& sequence, size_t start = 0,
			  size_t length = std::string::npos):
		_sequence(sequence), _start(start), _length(length)
	{}

	KmerIterator begin()
	{
		if (_sequence.length() < Parameters::get().kmerSize + _start)
			return this->end();

		return KmerIterator(&_sequence, _start);
	}

	KmerIterator end()
	{
		size_t end = _length == std::string::npos ?
						_sequence.length() : _length + _start;
		return KmerIterator(&_sequence, end - Parameters::get().kmerSize);
	}

private:
	const DnaSequence& _sequence;
	const size_t _start;
	const size_t _length;
};
