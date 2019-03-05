//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <cassert>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "flat_hash_map/flat_hash_map.hpp"

//Immutable dna sequence class
class DnaSequence
{
public:
	typedef size_t NuclType;

private:
	static const int NUCL_BITS = 2;
	static const int NUCL_IN_CHUNK = sizeof(NuclType) * 8 / NUCL_BITS;

	struct SharedBuffer
	{
		SharedBuffer(): useCount(0), length(0) {}
		size_t useCount;
		size_t length;
		std::vector<size_t> chunks;
	};

public:
	DnaSequence():
		_complement(false)
	{
		_data = new SharedBuffer;
		++_data->useCount;
	}

	~DnaSequence()
	{
		if (_data != nullptr)
		{
			//std::cout << "Destructor!\n";
			--_data->useCount;
			if (_data->useCount == 0) 
			{
				//std::cout << "Deleting\n";
				delete _data;
			}
		}
	}

	explicit DnaSequence(const std::string& string):
		_complement(false)
	{
		_data = new SharedBuffer;
		++_data->useCount;

		if (string.empty()) return;

		_data->length = string.length();
		_data->chunks.assign((_data->length - 1) / NUCL_IN_CHUNK + 1, 0);
		for (size_t i = 0; i < string.length(); ++i)
		{
			size_t chunkId = i / NUCL_IN_CHUNK;
			_data->chunks[chunkId] |= dnaToId(string[i]) << (i % NUCL_IN_CHUNK) * 2;
		}
	}

	DnaSequence(const DnaSequence& other):
		_data(other._data),
		_complement(other._complement)
	{
		++_data->useCount;
	}

	DnaSequence(DnaSequence&& other):
		_data(other._data),
		_complement(other._complement)
	{
		other._data = nullptr;
	}

	DnaSequence& operator=(const DnaSequence& other)
	{
		--_data->useCount;
		if (_data->useCount == 0) delete _data;

		_complement = other._complement;
		_data = other._data;
		++_data->useCount;
		return *this;
	}

	DnaSequence& operator=(DnaSequence&& other)
	{
		--_data->useCount;
		if (_data->useCount == 0) delete _data;

		_data = other._data;
		_complement = other._complement;
		other._data = nullptr;
		return *this;
	}

	size_t length() const {return _data->length;}

	char at(size_t index) const 
	{
		if (_complement)
		{
			index = _data->length - index - 1;
		}
		size_t id = (_data->chunks[index / NUCL_IN_CHUNK] >> 
					 (index % NUCL_IN_CHUNK) * 2 ) & 3;
		return idToDna(!_complement ? id : ~id & 3);
	}

	NuclType atRaw(size_t index) const 
	{
		if (_complement)
		{
			index = _data->length - index - 1;
		}
		size_t id = (_data->chunks[index / NUCL_IN_CHUNK] >> 
					 (index % NUCL_IN_CHUNK) * 2 ) & 3;
		return !_complement ? id : ~id & 3;
	}
	
	//TODO: use the same shared buffer
	
	DnaSequence complement() const
	{
		DnaSequence complSequence(*this);
		complSequence._complement = true;
		return complSequence;
	}

	DnaSequence substr(size_t start, size_t length) const;
	std::string str() const;	

	static size_t dnaToId(char c)
	{
		return _dnaTable[(size_t)c];
	}

	static char idToDna(size_t id)
	{
		static char table[] = {'A', 'C', 'G', 'T'};
		return table[id];
	}

private:
	static std::vector<size_t> _dnaTable;

	struct TableFiller
	{
		TableFiller()
		{
			static bool tableFilled = false;
			if (!tableFilled)
			{
				tableFilled = true;
				_dnaTable.assign(256, -1);	//256 chars
				_dnaTable[(size_t)'A'] = 0;
				_dnaTable[(size_t)'a'] = 0;
				_dnaTable[(size_t)'C'] = 1;
				_dnaTable[(size_t)'c'] = 1;
				_dnaTable[(size_t)'G'] = 2;
				_dnaTable[(size_t)'g'] = 2;
				_dnaTable[(size_t)'T'] = 3;
				_dnaTable[(size_t)'t'] = 3;
			}
		}
	};
	static TableFiller _filler;

	SharedBuffer* _data;
	bool _complement;
};

inline std::string DnaSequence::str() const 
{
	std::string result;
	result.reserve(this->length());
	for (size_t i = 0; i < this->length(); ++i)
	{
		result.push_back(this->at(i));
	}
	return result;
}

inline DnaSequence DnaSequence::substr(size_t start, size_t length) const 
{
	if (length == 0) throw std::runtime_error("Zero length subtring");
	if (start >= _data->length) throw std::runtime_error("Incorrect substring start");

	if (start + length > _data->length)
	{
		length = _data->length - start;
	}

	DnaSequence newSequence;
	newSequence._data->length = length;
	newSequence._data->chunks.assign((length - 1) / NUCL_IN_CHUNK + 1, 0);

	for (size_t i = 0; i < length; ++i)
	{
		size_t nucId = this->atRaw(start + i);
		size_t newChunkId = i / NUCL_IN_CHUNK;
		newSequence._data->chunks[newChunkId] |= nucId << (i % NUCL_IN_CHUNK) * 2;
	}

	return newSequence;
}
