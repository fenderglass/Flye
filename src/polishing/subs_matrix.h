//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>

typedef int64_t AlnScoreType;

class SubstitutionMatrix 
{
public:
	SubstitutionMatrix(const std::string& path);
	AlnScoreType getScore(char v, char w) const
	{
		return _matrix[(size_t)v * MAX_CHAR + (size_t)w];
	}
	
private:
	void loadMatrix(const std::string& path);
	void setScore(char v, char w, AlnScoreType score)
	{
		_matrix[(size_t)v * MAX_CHAR + (size_t)w] = score;
	}

	const size_t MAX_CHAR = std::numeric_limits<char>::max();
	std::vector<AlnScoreType> _matrix;
};

class HopoMatrix
{
public:
	//State represents a homopolymer that is observed in
	//reference sequence
	struct State
	{
		State():
			nucl(0), length(0), id(0)
		{}

		State(char nucl, uint32_t length);

		State(const std::string& str, size_t start = 0,
			  size_t end = std::string::npos);

		char nucl;
		uint32_t length;
		uint32_t id;
	};
	//Observation represents the read segment that corresponds
	//to a homopolymer in the reference (State). Might not be
	//a homopolymer, e.g. contain some other nucleotides, like
	//5A2X
	struct Observation
	{
		Observation(uint32_t id, bool extactMatch = false):
			id(id), extactMatch(extactMatch)
		{}
		uint32_t id;
		bool extactMatch;
	};
	typedef std::vector<Observation> ObsVector;

	HopoMatrix(const std::string& fileName);
	AlnScoreType getObsProb(State state, Observation observ) const
		{return _observationProbs[state.id][observ.id];}
	AlnScoreType getGenomeProb(State state) const
		{return _genomeProbs[state.id];}
	ObsVector knownObservations(State state) const;
	static Observation strToObs(char mainNucl, const std::string& dnaStr, 
								size_t start = 0, 
								size_t end = std::string::npos);

	//static std::string obsToStr(Observation obs);
private:
	void loadMatrix(const std::string& filaName);

	std::vector<std::vector<AlnScoreType>> _observationProbs;
	std::vector<AlnScoreType> 			   _genomeProbs;
};
