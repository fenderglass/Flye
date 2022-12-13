//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>

#include "../sequence/sequence_container.h"
#include "../common/config.h"
#include "../common/logger.h"
#include "../common/utils.h"
#include "../common/memory_info.h"

#include "../repeat_graph/repeat_graph.h"
#include "../repeat_graph/multiplicity_inferer.h"
#include "../repeat_graph/haplotype_resolver.h"
#include "../repeat_graph/graph_processing.h"
#include "../repeat_graph/repeat_resolver.h"
#include "../repeat_graph/output_generator.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outFolder, std::string& logFile, 
			   std::string& inAssembly, int& kmerSize,
			   int& minOverlap, bool& debug, size_t& numThreads, 
			   std::string& configPath, bool& unevenCov,
			   bool& keepHaplotypes, std::string& extraParams)
{
	auto printUsage = []()
	{
		std::cerr << "Usage: flye-repeat "
				  << " --disjointigs path --reads path --out-dir path --config path\n"
				  << "\t\t[--log path] [--treads num] [--kmer size] [--meta] [--keep-haplotypes]\n"
				  << "\t\t[--min-ovlp size] [--extra-params] [--debug] [-h]\n\n"
				  << "Required arguments:\n"
				  << "  --disjointigs path\tpath to disjointigs file\n"
				  << "  --reads path\tcomma-separated list of read files\n"
				  << "  --out-dir path\tpath to output directory\n"
				  << "  --config path\tpath to the config file\n\n"
				  << "Optional arguments:\n"
				  << "  --kmer size\tk-mer size [default = 15] \n"
				  << "  --min-ovlp size\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "  --debug \t\tenable debug output "
				  << "[default = false] \n"
				  << "  --meta \t\tenable uneven coverage (metagenome) mode "
				  << "[default = false] \n"
				  << "  --keep-haplotypes \t\tdo not collapse alternative haplotypes "
				  << "[default = false] \n"
				  << "  --log log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "  --extra-params additional config parameters "
				  << "[default = not set] \n"
				  << "  --threads num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};
	
	int optionIndex = 0;
	static option longOptions[] =
	{
		{"disjointigs", required_argument, 0, 0},
		{"reads", required_argument, 0, 0},
		{"out-dir", required_argument, 0, 0},
		{"config", required_argument, 0, 0},
		{"log", required_argument, 0, 0},
		{"threads", required_argument, 0, 0},
		{"kmer", required_argument, 0, 0},
		{"min-ovlp", required_argument, 0, 0},
		{"extra-params", required_argument, 0, 0},
		{"meta", no_argument, 0, 0},
		{"keep-haplotypes", no_argument, 0, 0},
		{"debug", no_argument, 0, 0},
		{0, 0, 0, 0}
	};

	int opt = 0;
	while ((opt = getopt_long(argc, argv, "h", longOptions, &optionIndex)) != -1)
	{
		switch(opt)
		{
		case 0:
			if (!strcmp(longOptions[optionIndex].name, "kmer"))
				kmerSize = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "threads"))
				numThreads = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "min-ovlp"))
				minOverlap = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "log"))
				logFile = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "debug"))
				debug = true;
			else if (!strcmp(longOptions[optionIndex].name, "meta"))
				unevenCov = true;
			else if (!strcmp(longOptions[optionIndex].name, "keep-haplotypes"))
				keepHaplotypes = true;
			else if (!strcmp(longOptions[optionIndex].name, "reads"))
				readsFasta = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "out-dir"))
				outFolder = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "disjointigs"))
				inAssembly = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "config"))
				configPath = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "extra-params"))
				extraParams = optarg;
			break;

		case 'h':
			printUsage();
			exit(0);
		}
	}
	if (readsFasta.empty() || outFolder.empty() || 
		inAssembly.empty() || configPath.empty())
	{
		printUsage();
		return false;
	}

	return true;
}

int repeat_main(int argc, char** argv)
{
	#ifdef NDEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	bool debugging = false;
	size_t numThreads = 1;
	int kmerSize = -1;
	int minOverlap = 5000;
	bool isMeta = false;
	bool keepHaplotypes = false; 
	std::string readsFasta;
	std::string inAssembly;
	std::string outFolder;
	std::string logFile;
	std::string configPath;
	std::string extraParams;
	if (!parseArgs(argc, argv, readsFasta, outFolder, logFile, inAssembly,
				   kmerSize, minOverlap, debugging, 
				   numThreads, configPath, isMeta, keepHaplotypes, extraParams))  return 1;
	
	Logger::get().setDebugging(debugging);
	if (!logFile.empty()) Logger::get().setOutputFile(logFile);
	Logger::get().debug() << "Build date: " << __DATE__ << " " << __TIME__;
	std::ios::sync_with_stdio(false);
	
	Logger::get().debug() << "Total RAM: " 
		<< getMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Available RAM: " 
		<< getFreeMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Total CPUs: " << std::thread::hardware_concurrency();

	Config::load(configPath);
	if (!extraParams.empty()) Config::addParameters(extraParams);
	if (kmerSize == -1)
	{
		kmerSize = Config::get("kmer_size");
	}
	Parameters::get().numThreads = numThreads;
	Parameters::get().kmerSize = kmerSize;
	Parameters::get().minimumOverlap = minOverlap;
	Parameters::get().unevenCoverage = isMeta;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 
	Logger::get().debug() << "Selected minimum overlap " << minOverlap;
	Logger::get().debug() << "Metagenome mode: " << "NY"[isMeta];

	Logger::get().info() << "Parsing disjointigs";
	SequenceContainer seqAssembly;
	std::vector<std::string> readsList = splitString(readsFasta, ',');
	try
	{
		seqAssembly.loadFromFile(inAssembly);
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	seqAssembly.buildPositionIndex();

	Logger::get().info() << "Building repeat graph";
	SequenceContainer edgeSequences;
	RepeatGraph rg(seqAssembly, &edgeSequences);
	rg.build();
	//rg.validateGraph();

	Logger::get().info() << "Parsing reads";
	SequenceContainer seqReads;
	try
	{
		for (auto& readsFile : readsList) seqReads.loadFromFile(readsFile);
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	seqReads.buildPositionIndex();
	//NOTE: it is important that we call the update below AFTER all reads
	//are loaded. It ensures that all sequences for the repeat
	//graph are stored in a continious chunk of memory.
	rg.updateEdgeSequences();

	Logger::get().info() << "Aligning reads to the graph";
	ReadAligner aligner(rg, seqReads);
	aligner.alignReads();
	MultiplicityInferer multInf(rg, aligner, seqAssembly);
	multInf.estimateCoverage();
	//aligner.storeAlignments(outFolder + "/read_alignment_before_rr");

	Logger::get().info() << "Simplifying the graph";

	multInf.removeUnsupportedEdges(/*only tips*/ true);
	//multInf.removeUnsupportedConnections();
	//rg.validateGraph();
	
	RepeatResolver repResolver(rg, seqAssembly, seqReads, aligner, multInf);
	HaplotypeResolver hapResolver(rg, aligner, seqAssembly, seqReads);
	GraphProcessor proc(rg, seqAssembly);
	OutputGenerator outGen(rg, aligner);

	//revealing/masking haplotypes
	auto markHaplotypes = [&hapResolver, isMeta]()
	{
		hapResolver.resetEdges();
		hapResolver.findHeterozygousLoops();
		hapResolver.findHeterozygousBulges();
		if (isMeta)
		{
			hapResolver.findRoundabouts();
			hapResolver.findSuperbubbles();
		}
	};

	//dump graph before first repeat resolution iteration
	markHaplotypes();
	repResolver.findRepeats();
	outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_before_rr.gv");
	outGen.outputFasta(proc.getEdgesPaths(), outFolder + "/graph_before_rr.fasta");
	if ((bool)Config::get("output_gfa_before_rr"))
	{
		outGen.outputGfa(proc.getEdgesPaths(), outFolder + "/graph_before_rr.gfa");
	}

	if (isMeta && !keepHaplotypes) 
	{
		repResolver.resolveSimpleRepeats();
	}
	for (int iterNum = 1; ;++iterNum)
	{
		int actions = 0;
		Logger::get().debug() << "[SIMPL] == Iteration " << iterNum << " ==";

		actions += multInf.splitNodes();
		if (isMeta) 
		{
			actions += multInf.disconnectMinorPaths();
		}
		actions += multInf.trimTips();

		//resolving repeats
		markHaplotypes();
		repResolver.findRepeats();
		actions += repResolver.resolveRepeats();

		//rg.validateGraph();

		if (!actions) break;
	}

	if (isMeta) 
	{
		multInf.resolveForks();
	}

	if (!keepHaplotypes)
	{
		hapResolver.collapseHaplotypes();
		repResolver.resolveSimpleRepeats();
	}
	else
	{
		//if we keeping haplotypes, update them one more time to ensure it's consistent
		//with the final graph structure
		markHaplotypes();
	}

	multInf.removeUnsupportedEdges(/*only tips*/ true);

	repResolver.findRepeats();
	repResolver.finalizeGraph();
	//rg.validateGraph();

	outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_after_rr.gv");
	rg.storeGraph(outFolder + "/repeat_graph_dump");
	aligner.storeAlignments(outFolder + "/read_alignment_dump");
	SequenceContainer::writeFasta(edgeSequences.iterSeqs(), 
								  outFolder + "/repeat_graph_edges.fasta",
								  /*only pos strand*/ true);

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	return 0;
}
