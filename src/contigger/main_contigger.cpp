//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <execinfo.h>

#include "../sequence/sequence_container.h"
#include "../common/config.h"
#include "../common/logger.h"
#include "../common/utils.h"
#include "../common/memory_info.h"

#include "../repeat_graph/repeat_graph.h"
#include "../repeat_graph/read_aligner.h"
#include "../repeat_graph/output_generator.h"
#include "../contigger/contig_extender.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outFolder, std::string& logFile, 
			   std::string& inGraphEdges, int& kmerSize,
			   int& minOverlap, bool& debug, size_t& numThreads, 
			   std::string& configPath, std::string& inRepeatGraph,
			   std::string& inReadsAlignment, bool& noScaffold,
			   std::string& extraParams)
{
	auto printUsage = []()
	{
		std::cerr << "Usage: flye-contigger "
				  << " --graph-edges path --reads path --out-dir path --config path\n"
				  << "\t\t--repeat-graph path --graph-aln path\n"
				  << "\t\t[--log path] [--treads num] [--kmer size] [--no-scaffold]\n"
				  << "\t\t[--min-ovlp size] [--debug] [--extra-params] [-h]\n\n"
				  << "Required arguments:\n"
				  << "  --graph-edges path\tpath to fasta with graph edges\n"
				  << "  --repeat-graph path\tpath to serialized repeat graph\n"
				  << "  --graph-aln path\tpath to read-graph alignment\n"
				  << "  --reads path\tcomma-separated list of read files\n"
				  << "  --out-dir path\tpath to output directory\n"
				  << "  --config path\tpath to the config file\n\n"
				  << "Optional arguments:\n"
				  << "  --kmer size\tk-mer size [default = 15] \n"
				  << "  --min-ovlp size\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "  --debug \t\tenable debug output "
				  << "[default = false] \n"
				  << "  --no-scaffold \t\tdisable scaffolding "
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
		{"graph-edges", required_argument, 0, 0},
		{"reads", required_argument, 0, 0},
		{"out-dir", required_argument, 0, 0},
		{"config", required_argument, 0, 0},
		{"repeat-graph", required_argument, 0, 0},
		{"graph-aln", required_argument, 0, 0},
		{"log", required_argument, 0, 0},
		{"threads", required_argument, 0, 0},
		{"kmer", required_argument, 0, 0},
		{"min-ovlp", required_argument, 0, 0},
		{"extra-params", required_argument, 0, 0},
		{"debug", no_argument, 0, 0},
		{"no-scaffold", no_argument, 0, 0},
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
			else if (!strcmp(longOptions[optionIndex].name, "no-scaffold"))
				noScaffold = true;
			else if (!strcmp(longOptions[optionIndex].name, "reads"))
				readsFasta = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "out-dir"))
				outFolder = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "graph-edges"))
				inGraphEdges = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "repeat-graph"))
				inRepeatGraph = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "graph-aln"))
				inReadsAlignment = optarg;
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
		inGraphEdges.empty() || configPath.empty() ||
		inRepeatGraph.empty() || inReadsAlignment.empty())
	{
		printUsage();
		return false;
	}

	return true;
}

int contigger_main(int argc, char** argv)
{
	#ifdef NDEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	bool debugging = false;
	size_t numThreads = 1;
	int kmerSize = -1;
	int minOverlap = 5000;
	bool noScaffold = false;
	std::string readsFasta;
	std::string inGraphEdges;
	std::string inRepeatGraph;
	std::string inReadsAlignment;
	std::string outFolder;
	std::string logFile;
	std::string configPath;
	std::string extraParams;
	if (!parseArgs(argc, argv, readsFasta, outFolder, logFile, inGraphEdges,
				   kmerSize, minOverlap, debugging, 
				   numThreads, configPath, inRepeatGraph, 
				   inReadsAlignment, noScaffold, extraParams))  return 1;
	
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
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 
	Logger::get().debug() << "Selected minimum overlap " << minOverlap;

	Logger::get().info() << "Reading sequences";
	SequenceContainer seqGraphEdges; 
	SequenceContainer seqReads;
	std::vector<std::string> readsList;
	readsList.push_back(readsFasta);
	try
	{
		seqGraphEdges.loadFromFile(inGraphEdges);
		for (auto& readsFile : readsList)
		{
			seqReads.loadFromFile(readsFile);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	seqReads.buildPositionIndex();
	//seqAssembly.buildPositionIndex();

	SequenceContainer emptyContainer;
	RepeatGraph rg(emptyContainer, &seqGraphEdges);
	rg.loadGraph(inRepeatGraph);
	//rg.validateGraph();
	ReadAligner aln(rg, seqReads);
	aln.loadAlignments(inReadsAlignment);
	OutputGenerator outGen(rg, aln);

	//Logger::get().info() << "Generating contigs";

	ContigExtender extender(rg, aln, emptyContainer, seqReads);
	extender.generateUnbranchingPaths();
	extender.generateContigs();
	extender.outputContigs(outFolder + "/contigs.fasta");
	extender.outputStatsTable(outFolder + "/contigs_stats.txt");

	std::string scaffoldFile = outFolder + "/scaffolds_links.txt"; 
	if (!noScaffold)
	{
		extender.outputScaffoldConnections(scaffoldFile);
	}
	else
	{
		std::ofstream scfFile(scaffoldFile);	//creates empty file
	}

	//outGen.dumpRepeats(extender.getUnbranchingPaths(),
	//				   outFolder + "/repeats_dump.txt");
	outGen.outputDot(extender.getUnbranchingPaths(),
					 outFolder + "/graph_final.gv");
	outGen.outputFasta(extender.getUnbranchingPaths(),
					   outFolder + "/graph_final.fasta");
	outGen.outputGfa(extender.getUnbranchingPaths(),
					 outFolder + "/graph_final.gfa");

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	return 0;
}
