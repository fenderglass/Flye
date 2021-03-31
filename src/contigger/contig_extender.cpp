//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "contig_extender.h"
#include "../repeat_graph/output_generator.h"
#include <cmath>

void ContigExtender::generateUnbranchingPaths()
{
	GraphProcessor proc(_graph, _asmSeqs);
	_unbranchingPaths = proc.getUnbranchingPaths();

	//tryna synchornize orientations of alternative
	//contigs from the same group
	std::unordered_set<FastaRecord::Id> flipped;
	for (auto& path : _unbranchingPaths)
	{
		int numEven = 0;
		for (auto& edge : path.path)
		{
			if (edge->altGroupId != -1 && edge->altGroupId % 2 == 0) ++numEven;
		}
		//trying to make groups with even alt ids to have positive path ids
		if (!path.id.strand() && numEven > (int)path.path.size() / 2)
		{
			flipped.insert(path.id);
			flipped.insert(path.id.rc());
		}
	}
	for (auto& path : _unbranchingPaths)
	{
		if (flipped.count(path.id))
		{
			path.id = path.id.rc();
		}
	}
	Logger::get().debug() << "Flipped " << flipped.size() / 2;

	_edgeToPath.clear();
	for (auto& path : _unbranchingPaths)
	{
		path.prefix = "edge_";
		if (path.id.strand())
		{
			Logger::get().debug() << "UPath " << path.id.signedId() 
							<< ": " << path.edgesStr();
		}

		for (auto& edge : path.path)
		{
			_edgeToPath[edge] = &path;
		}
	}

	Logger::get().debug() << "Final graph contain " 
		<< _unbranchingPaths.size() / 2 << " egdes";
}


void ContigExtender::generateContigs()
{
	Logger::get().debug() << "Extending contigs into repeats";

	bool graphContinue = (bool)Config::get("extend_contigs_with_repeats");

	OutputGenerator outGen(_graph, _aligner);
	auto coreSeqs = outGen.generatePathSequences(_unbranchingPaths);
	std::unordered_map<UnbranchingPath*, FastaRecord*> upathsSeqs;
	for (size_t i = 0; i < _unbranchingPaths.size(); ++i)
	{
		upathsSeqs[&_unbranchingPaths[i]] = &coreSeqs[i];
	}

	std::unordered_map<GraphEdge*, 
					   std::vector<const GraphAlignment*>> alnIndex;
	for (auto& aln : _aligner.getAlignments())
	{
		if (aln.size() > 1)
		{
			for (auto edge : aln)
			{
				alnIndex[edge.edge].push_back(&aln);
			}
		}
	}

	std::unordered_set<GraphEdge*> coveredRepeats;
	std::unordered_map<const GraphEdge*, bool> repeatDirections;
	auto canTraverse = [&repeatDirections] (const GraphEdge* edge)
	{
		//if (edge->isLooped() && edge->selfComplement) return false;
		return !repeatDirections.count(edge) || 
			   repeatDirections.at(edge);
	};

	typedef std::pair<GraphPath, std::string> PathAndSeq;
	auto extendPathRight =
		[this, &coveredRepeats, &repeatDirections, &upathsSeqs, 
		 &canTraverse, &alnIndex, graphContinue] 
	(UnbranchingPath& upath)
	{

		bool extendFwd = !upath.path.back()->nodeRight->outEdges.empty();
		if (!extendFwd) return PathAndSeq();

		//first, choose the longest aligned read from this edge
		int32_t maxExtension = 0;
		GraphAlignment bestAlignment;
		for (auto pathPtr : alnIndex[upath.path.back()])
		{
			const GraphAlignment& path = *pathPtr;
			for (size_t i = 0; i < path.size(); ++i)
			{
				if (path[i].edge == upath.path.back() &&
					i < path.size() - 1)
				{
					size_t j = i + 1;
					while (j < path.size() && 
						   path[j].edge->repetitive &&
						   !path[j].edge->altHaplotype &&
						   canTraverse(path[j].edge)) ++j;
					if (j == i + 1) break;

					int32_t alnLen = path[j - 1].overlap.curEnd - 
									 path[i + 1].overlap.curBegin;
					if (alnLen > maxExtension)
					{
						maxExtension = alnLen;
						bestAlignment.clear();
						std::copy(path.begin() + i + 1, path.begin() + j,
								  std::back_inserter(bestAlignment));
					}
					break;
				}
			}
		}
		if (maxExtension == 0) return PathAndSeq();

		auto upathAln = this->asUpathAlignment(bestAlignment);
		auto lastUpath = upathAln.back().upath;
		int32_t overhang = upathsSeqs[lastUpath]->sequence.length() - 
						   upathAln.back().aln.back().overlap.curEnd + 
						   upathAln.back().aln.front().overlap.curBegin;
		bool lastIncomplete = overhang > (int)Config::get("max_separation");
		//Logger::get().debug() << "Ctg " << upath.id.signedId() <<
		//	" overhang " << overhang << " upath " << lastUpath->id.signedId();

		for (size_t i = 0; i < upathAln.size(); ++i)
		{
			//in case we don't extend using graph structure,
			//dont mark the last upath as covered if it is
			//incomplete
			if (i == upathAln.size() - 1 && 
				lastIncomplete && !graphContinue) break;

			for (auto& aln : upathAln[i].aln)
			{
				repeatDirections[aln.edge] = true;
				repeatDirections[_graph.complementEdge(aln.edge)] = false;
				coveredRepeats.insert(aln.edge);
				coveredRepeats.insert(_graph.complementEdge(aln.edge));
			}
		}

		//generate extension sequence
		std::string extendedSeq;
		if (lastIncomplete && graphContinue)
		{
			upathAln.pop_back();
		}
		if (!upathAln.empty())
		{
			FastaRecord::Id readId = bestAlignment.front().overlap.curId;
			int32_t readStart = upathAln.front().aln.front().overlap.curBegin;
			int32_t readEnd = upathAln.back().aln.back().overlap.curEnd;
			extendedSeq = _readSeqs.getSeq(readId)
				.substr(readStart, readEnd - readStart).str();
		}
		if (lastIncomplete && graphContinue)
		{
			extendedSeq += upathsSeqs[lastUpath]->sequence.str();
		}
		
		GraphPath extendedPath;
		for (auto& ualn : upathAln)
		{
			for (auto& edgeAln : ualn.aln) extendedPath.push_back(edgeAln.edge);
		}
		if (lastIncomplete && graphContinue)
		{
			for (auto& edge : lastUpath->path) extendedPath.push_back(edge);
		}
		return PathAndSeq(extendedPath, extendedSeq);
	};

	std::unordered_map<FastaRecord::Id, UnbranchingPath*> idToPath;
	for (auto& ctg : _unbranchingPaths)
	{
		idToPath[ctg.id] = &ctg;
	}

	for (auto& upath : _unbranchingPaths)
	{
		if (upath.repetitive || !upath.id.strand()) continue;
		if (!idToPath.count(upath.id.rc())) continue;	//self-complement

		auto rightExt = extendPathRight(upath);
		auto leftExt = extendPathRight(*idToPath[upath.id.rc()]);
		leftExt.first = _graph.complementPath(leftExt.first);
		leftExt.second = DnaSequence(leftExt.second).complement().str();

		Contig contig(upath);
		auto leftPaths = this->asUpaths(leftExt.first);
		auto rightPaths = this->asUpaths(rightExt.first);

		GraphPath leftEdges;
		for (auto& path : leftPaths)
		{
			leftEdges.insert(leftEdges.end(), path->path.begin(), 
						     path->path.end());
		}
		contig.graphEdges.path.insert(contig.graphEdges.path.begin(), 
								      leftEdges.begin(), leftEdges.end());
		contig.graphPaths.insert(contig.graphPaths.begin(), 
								 leftPaths.begin(), leftPaths.end());

		GraphPath rightEdges;
		for (auto& path : rightPaths)
		{
			rightEdges.insert(rightEdges.end(), path->path.begin(), 
						      path->path.end());
		}
		contig.graphEdges.path.insert(contig.graphEdges.path.end(), 
								      rightEdges.begin(), rightEdges.end());
		contig.graphPaths.insert(contig.graphPaths.end(), 
								 rightPaths.begin(), rightPaths.end());

		auto coreSeq = upathsSeqs[&upath]->sequence.str();
		contig.sequence = DnaSequence(leftExt.second + coreSeq + rightExt.second);

		_contigs.push_back(std::move(contig));
	}

	//add repetitive contigs that were not covered by the extended paths
	int numCovered = 0;
	for (auto& upath : _unbranchingPaths)
	{
		if (!upath.repetitive || !upath.id.strand()) continue;

		bool covered = false;
		for (auto& edge : upath.path)
		{
			if (coveredRepeats.count(edge)) covered = true;
		}
		if (!covered)
		{
			_contigs.emplace_back(upath);
			_contigs.back().sequence = upathsSeqs[&upath]->sequence;
		}
		else
		{
			++numCovered;
			//Logger::get().debug() << "Covered: " << upath.id.signedId();
		}
	}

	for (auto& ctg : _contigs)
	{
		ctg.graphEdges.prefix = "contig_";
	}

	Logger::get().debug() << "Covered " << numCovered << " repetitive contigs";
	Logger::get().info() << "Generated " << _contigs.size() << " contigs";
}

std::vector<UnbranchingPath*> ContigExtender::asUpaths(const GraphPath& path)
{
	std::vector<UnbranchingPath*> upathRepr;
	for (size_t i = 0; i < path.size(); ++i)
	{
		UnbranchingPath* upath = _edgeToPath.at(path[i]);
		if (upathRepr.empty() || upathRepr.back() != upath ||
			path[i - 1] == path[i])
		{
			upathRepr.push_back(upath);
		}
	}

	return upathRepr;
}

std::vector<ContigExtender::UpathAlignment> 
	ContigExtender::asUpathAlignment(const GraphAlignment& graphAln)
{
	std::vector<ContigExtender::UpathAlignment> upathAln;
	for (size_t i = 0; i < graphAln.size(); ++i)
	{
		UnbranchingPath* upath = _edgeToPath.at(graphAln[i].edge);
		if (upathAln.empty() || upathAln.back().upath != upath ||
			graphAln[i - 1].edge == graphAln[i].edge)
		{
			upathAln.emplace_back();
			upathAln.back().upath = upath;
		}
		upathAln.back().aln.push_back(graphAln[i]);
	}

	return upathAln;
}

void ContigExtender::outputStatsTable(const std::string& filename)
{
	std::ofstream fout(filename);
	if (!fout) throw std::runtime_error("Can't write " + filename);

	fout << "#seq_name\tlength\tcoverage\tcircular\trepeat"
		<< "\tmult\ttelomere\talt_group\tgraph_path\n";

	char YES_NO[] = {'N', 'Y'};

	//TODO: compute mean coverage
	int64_t sumCov = 0;
	int64_t sumLength = 0;
	for (auto& edge : _graph.iterEdges())
	{
		if (edge->edgeId.strand())
		{
			sumCov += edge->meanCoverage * edge->length();
			sumLength += edge->length();
		}
	}
	int meanCoverage = sumCov / (sumLength + 1);

	for (auto& ctg : _contigs)
	{
		std::string pathStr;
		for (auto& upath : ctg.graphPaths)
		{
			pathStr += std::to_string(upath->id.signedId()) + ",";
		}
		pathStr.pop_back();

		int estMult = std::max(1.0f, std::round((float)ctg.graphEdges.meanCoverage / 
											    meanCoverage));

		int altGroup = ctg.graphEdges.path.front()->altGroupId;
		for (auto& edge : ctg.graphEdges.path)
		{
			if (edge->altGroupId != altGroup)
			{
				altGroup = -1;
				break;
			}
		}
		std::string altGroupStr = (altGroup == -1) ? "*" : std::to_string(altGroup);

		std::string telomereStr;
		bool telLeft = (ctg.graphEdges.path.front()->nodeLeft->isTelomere());
		bool telRight = (ctg.graphEdges.path.back()->nodeRight->isTelomere());
		if (telLeft && !telRight) telomereStr = "left";
		if (!telLeft && telRight) telomereStr = "right";
		if (telLeft && telRight) telomereStr = "both";
		if (!telLeft && !telRight) telomereStr = "none";

		fout << ctg.graphEdges.name() << "\t" << ctg.sequence.length() << "\t" 
			<< ctg.graphEdges.meanCoverage << "\t"
			<< YES_NO[ctg.graphEdges.circular]
			<< "\t" << YES_NO[ctg.graphEdges.repetitive] << "\t"
			<< estMult << "\t" << telomereStr << "\t" << altGroupStr << "\t"
			<< pathStr << "\n";

		//Logger::get().debug() << "Contig: " << ctg.graphEdges.id.signedId()
		//	<< ": " << pathStr;
	}
}

void ContigExtender::outputContigs(const std::string& filename)
{
	std::vector<FastaRecord> contigsFasta;
	for (auto& ctg : _contigs)
	{
		contigsFasta.emplace_back(ctg.sequence, ctg.graphEdges.name(), 
					   			  FastaRecord::ID_NONE);
	}
	SequenceContainer::writeFasta(contigsFasta, filename);
}

void ContigExtender::outputScaffoldConnections(const std::string& filename)
{
	Logger::get().debug() << "Generating scaffold connections";

	std::ofstream fout(filename);
	if (!fout) throw std::runtime_error("Can't open " + filename);

	auto reachableEdges = [this](GraphEdge* startEdge)
	{
		std::vector<GraphEdge*> dfsStack;
		std::unordered_set<GraphEdge*> visited;
		std::unordered_set<GraphEdge*> reachableUnique;

		dfsStack.push_back(startEdge);
		while(!dfsStack.empty())
		{
			GraphEdge* curEdge = dfsStack.back(); 
			dfsStack.pop_back();
			if (visited.count(curEdge)) continue;

			visited.insert(curEdge);
			visited.insert(_graph.complementEdge(curEdge));
			for (auto& adjEdge: curEdge->nodeRight->outEdges)
			{
				if (adjEdge->isRepetitive() && !adjEdge->isLooped() &&
					!visited.count(adjEdge))
				{
					dfsStack.push_back(adjEdge);
				}
				else if (!adjEdge->isRepetitive() &&
						 adjEdge->edgeId != startEdge->edgeId &&
						 adjEdge->edgeId != startEdge->edgeId.rc())
				{
					reachableUnique.insert(adjEdge);
					if (reachableUnique.size() > 1) return reachableUnique;
				}
			}
		}
		return reachableUnique;
	};

	int numScaffolds = 0;
	for (auto& edge : _graph.iterEdges())
	{
		if (edge->repetitive) continue;

		auto reachableUnique = reachableEdges(edge);	
		if (reachableUnique.size() != 1) continue;
		GraphEdge* outEdge = *reachableUnique.begin();
		if (reachableEdges(_graph.complementEdge(outEdge)).size() != 1) continue;

		if ((edge->nodeRight->isBifurcation() || 
			 outEdge->nodeLeft->isBifurcation()) &&
			edge->edgeId != outEdge->edgeId.rc() &&
			abs(edge->edgeId.signedId()) < abs(outEdge->edgeId.signedId()))
		{
			UnbranchingPath leftCtg = *this->asUpaths({edge}).front();
			UnbranchingPath rightCtg = *this->asUpaths({outEdge}).front();
			if (leftCtg.id != rightCtg.id)
			{
				leftCtg.prefix = "contig_";
				rightCtg.prefix = "contig_";
				fout << leftCtg.nameUnsigned() << "\t" << 
					(leftCtg.id.strand() ? '+' : '-') << "\t" <<
					rightCtg.nameUnsigned() << "\t" << 
					(rightCtg.id.strand() ? '+' : '-') << "\n";
				++numScaffolds;
			}
		}
	}
	Logger::get().info() << "Added " << numScaffolds << " scaffold connections";
}
