//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <algorithm>

#include "homo_polisher.h"
#include "../common/matrix.h"
#include "flat_hash_map/flat_hash_map.hpp"


namespace
{
	//Computes global pairwise alignment with custom substitution matrix
	void pairwiseAlignment(const std::string& seqOne, const std::string& seqTwo,
						   const SubstitutionMatrix& subsMat,
						   std::string& outOne, std::string& outTwo)
	{
		Matrix<AlnScoreType> scoreMat(seqOne.length() + 1, seqTwo.length() + 1);
		Matrix<char> backtrackMat(seqOne.length() + 1, seqTwo.length() + 1);

		scoreMat.at(0, 0) = 0;
		backtrackMat.at(0, 0) = 0;
		for (size_t i = 0; i < seqOne.length(); ++i) 
		{
			scoreMat.at(i + 1, 0) = scoreMat.at(i, 0) + 
									subsMat.getScore(seqOne[i], '-');
			backtrackMat.at(i + 1, 0) = 1;
		}
		for (size_t i = 0; i < seqTwo.length(); ++i) 
		{
			scoreMat.at(0, i + 1) = scoreMat.at(0, i) + 
									subsMat.getScore('-', seqTwo[i]);
			backtrackMat.at(0, i + 1) = 0;
		}

		//filling DP matrices
		for (size_t i = 1; i < seqOne.length() + 1; ++i)
		{
			for (size_t j = 1; j < seqTwo.length() + 1; ++j) 
			{
				AlnScoreType left = scoreMat.at(i, j - 1) + 
							 		subsMat.getScore('-', seqTwo[j - 1]);
				AlnScoreType up = scoreMat.at(i - 1, j) + 
								  subsMat.getScore(seqOne[i - 1], '-');
				AlnScoreType cross = scoreMat.at(i - 1, j - 1) + 
							  		 subsMat.getScore(seqOne[i - 1], seqTwo[j - 1]);

				int prev = 2;
				AlnScoreType score = cross;
				if (up > score)
				{
					prev = 1;
					score = up;
				}
				if (left > score)
				{
					prev = 0;
					score = left;
				}
				scoreMat.at(i, j) = score;
				backtrackMat.at(i, j) = prev;
			}
		}

		//backtrack
		int i = seqOne.length();
		int j = seqTwo.length();
		outOne.clear();
		outTwo.clear();

		while (i != 0 || j != 0) 
		{
			if(backtrackMat.at(i, j) == 1) 
			{
				outOne += seqOne[i - 1];
				outTwo += '-';
				i -= 1;
			}
			else if (backtrackMat.at(i, j) == 0) 
			{
				outOne += '-';
				outTwo += seqTwo[j - 1];
				j -= 1;
			}
			else
			{
				outOne += seqOne[i - 1];
				outTwo += seqTwo[j - 1];
				i -= 1;
				j -= 1;
			}
		}
		std::reverse(outOne.begin(), outOne.end());
		std::reverse(outTwo.begin(), outTwo.end());
		outOne += "$";
		outTwo += "$";
	}

	//Splits aligned strings into homopolymer runs (wrt to candAln)
	std::vector<std::pair<HopoMatrix::State, HopoMatrix::Observation>>
	splitBranchHopos(const std::string& candAln, const std::string& branchAln)
	{
		//std::cerr << candAln << std::endl << branchAln << std::endl << std::endl;
		std::vector<std::pair<HopoMatrix::State, 
							  HopoMatrix::Observation>> result;
		size_t prevPos = 0;
		while (candAln[prevPos] == '-') ++prevPos;
		char prevNucl = candAln[prevPos];

		int runLength = 0;
		int gapLength = 0;
		int hopoCount = 0;
		for (size_t pos = prevPos + 1; pos < candAln.length(); ++pos)
		{
			if (candAln[pos] != '-') ++runLength;
			if (candAln[pos] == '-')
			{
				++gapLength;
			}
			else 
			{
				if (candAln[pos] != prevNucl)
				{
					bool leftMatch = prevPos == 0 || 
							candAln[prevPos - 1] == branchAln[prevPos - 1] ||
							branchAln[prevPos - 1] == candAln[prevPos];
					bool rightMatch = pos == candAln.length() - 1 || 
							candAln[pos] == branchAln[pos] ||
							branchAln[pos] == candAln[pos - 1];

					++hopoCount;
					//wobble
					size_t branchPrevPos = prevPos -
						(size_t)(prevPos > 0 && branchAln[prevPos - 1] == 
												candAln[prevPos]);
					size_t branchPos = pos + 
						(size_t)(branchAln[pos] == candAln[pos - 1]);

					/*if (prevPos > 0 && prevPos < branchAln.length() - 1)
					{
						std::string branchSubseq(branchAln, prevPos - 1, pos - prevPos + 2);
						std::string candSubseq(candAln, prevPos - 1, pos - prevPos + 2);
						std::cout << candSubseq << " " << std::endl 
								  << branchSubseq << std::endl 
								  << (rightMatch && leftMatch) << std::endl;
					}*/

					auto state = HopoMatrix::State(candAln, prevPos, pos);
					auto observ = HopoMatrix::strToObs(state.nucl, branchAln, 
													   branchPrevPos, branchPos);
					observ.extactMatch = leftMatch && rightMatch;
					result.emplace_back(state, observ);

					prevNucl = candAln[pos];
					prevPos = pos - gapLength;
					gapLength = 0;
					runLength = 0;
				}
				else
				{
					gapLength = 0;
				}
			}
		}

		return result;
	}
}

//processes a single bubble
void HomoPolisher::polishBubble(Bubble& bubble) const
{
	std::string prevCandidate;
	std::string curCandidate = bubble.candidate;

	std::vector<HopoMatrix::State> states;
	std::vector<HopoMatrix::ObsVector> observations;

	for (auto& branch : bubble.branches)
	{
		std::string alnCand;
		std::string alnBranch;
		pairwiseAlignment(bubble.candidate, branch, _subsMatrix,
						  alnCand, alnBranch);

		auto splitHopo = splitBranchHopos(alnCand, alnBranch);
		if (states.empty())
		{
			states.assign(splitHopo.size(), HopoMatrix::State());
			observations.assign(splitHopo.size(), HopoMatrix::ObsVector());
		}
		assert(states.size() == splitHopo.size());

		for (size_t i = 0; i < splitHopo.size(); ++i)
		{
			states[i] = splitHopo[i].first;
			observations[i].push_back(splitHopo[i].second);
		}
	}

	std::string newConsensus;
	for (size_t i = 0; i < states.size(); ++i)
	{
		size_t length = states[i].length;
		if (length > 1)	//only homopolymers
		{
			length = this->mostLikelyLen(states[i].nucl, observations[i]);
		}
		newConsensus += std::string(length, states[i].nucl);
		/*if (length != (size_t)states[i].length)
		{

			std::cout << (int)states[i].length << states[i].nucl 
					  << " -> " << length << std::endl;
		}*/
	}

	if (newConsensus != bubble.candidate)
	{
		StepInfo info;
		info.sequence = newConsensus;
		bubble.polishSteps.push_back(info);
		bubble.candidate = newConsensus;
	}
}

//likelihood of a given state
AlnScoreType HomoPolisher::likelihood(HopoMatrix::State state, 
								const HopoMatrix::ObsVector& observations) const
{
	AlnScoreType likelihood = 0.0f;
	for (auto obs : observations)
	{
		if (obs.extactMatch)
		{
			likelihood += _hopoMatrix.getObsProb(state, obs);
		}
	}
	likelihood += _hopoMatrix.getGenomeProb(state);
	return likelihood;
}

//for a given homopolymer and set of observations,
//computes most likely homopolymer length
size_t HomoPolisher::mostLikelyLen(char nucleotide,
								   const HopoMatrix::ObsVector& 
								   observations) const
{
	assert(!observations.empty());
	const size_t MIN_HOPO = 1;
	const size_t MAX_HOPO = 20;

	typedef std::pair<AlnScoreType, size_t> ScorePair;
	std::vector<ScorePair> scores;
	for (size_t len = MIN_HOPO; len <= MAX_HOPO; ++len)
	{
		auto newState = HopoMatrix::State(nucleotide, len);
		AlnScoreType likelihood = this->likelihood(newState, observations);
		scores.push_back(std::make_pair(likelihood, len));
	}

	std::sort(scores.begin(), scores.end(), 
			  [](const ScorePair& p1, const ScorePair& p2)
			  {return p1.first > p2.first;});

	size_t maxRun = this->compareTopTwo(nucleotide, scores[0].second, 
										scores[1].second, observations);
	return maxRun;
}

//Compares top two homopolimer candidates in a more precise manner
size_t HomoPolisher::compareTopTwo(char nucleotide, size_t firstChoice, 
								   size_t secondChoice,
				   				   const HopoMatrix::ObsVector& 
								   observations) const
{
	size_t choices[] = {firstChoice, secondChoice};
	HopoMatrix::ObsVector knownObs[2];

	for (size_t i = 0; i < 2; ++i)
	{
		auto state = HopoMatrix::State(nucleotide, choices[i]);
		knownObs[i] = _hopoMatrix.knownObservations(state);
	}
	
	//getting common known observations
	ska::flat_hash_set<uint32_t> fstSet;
	for (auto obs : knownObs[0]) fstSet.insert(obs.id);
	ska::flat_hash_set<uint32_t> commonSet;
	for (auto obs : knownObs[1])
	{
		if (fstSet.count(obs.id)) commonSet.insert(obs.id);
	}
	HopoMatrix::ObsVector commonObservations;
	for (auto obs : observations)
	{
		if (commonSet.count(obs.id)) commonObservations.push_back(obs);
	}

	AlnScoreType likelihoods[2];
	for (size_t i = 0; i < 2; ++i)
	{
		auto state = HopoMatrix::State(nucleotide, choices[i]);
		likelihoods[i] = this->likelihood(state, commonObservations);
	}

	return (likelihoods[0] > likelihoods[1]) ? choices[0] : choices[1];
}
