//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <execinfo.h>

#include "../sequence/vertex_index.h"
#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../sequence/consensus_generator.h"
#include "../common/config.h"
#include "../assemble/extender.h"
#include "../assemble/parameters_estimator.h"
#include "../common/logger.h"
#include "../common/utils.h"
#include "../common/memory_info.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outAssembly, std::string& logFile, size_t& genomeSize,
			   int& kmerSize, bool& debug, size_t& numThreads, int& minOverlap, 
			   std::string& configPath, int& minReadLength, bool& unevenCov, 
			   std::string& extraParams)
{
	auto printUsage = []()
	{
		std::cerr << "Usage: flye-assemble "
				  << " --reads path --out-asm path --config path [--genome-size size]\n"
				  << "\t\t[--min-read length] [--log path] [--treads num] [--extra-params]\n"
				  << "\t\t[--kmer size] [--meta] [--min-ovlp size] [--debug] [-h]\n\n"
				  << "Required arguments:\n"
				  << "  --reads path\tcomma-separated list of read files\n"
				  << "  --out-asm path\tpath to output file\n"
				  << "  --config path\tpath to the config file\n\n"
				  << "Optional arguments:\n"
				  << "  --genome-size size\tgenome size in bytes\n"
				  << "  --kmer size\tk-mer size [default = 15] \n"
				  << "  --min-ovlp size\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "  --debug \t\tenable debug output "
				  << "[default = false] \n"
				  << "  --meta \t\tenable uneven coverage (metagenome) mode "
				  << "[default = false] \n"
				  << "  --extra-params additional config parameters "
				  << "[default = not set] \n"
				  << "  --log log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "  --threads num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};
	
	int optionIndex = 0;
	static option longOptions[] =
	{
		{"reads", required_argument, 0, 0},
		{"out-asm", required_argument, 0, 0},
		{"genome-size", required_argument, 0, 0},
		{"config", required_argument, 0, 0},
		{"min-read", required_argument, 0, 0},
		{"log", required_argument, 0, 0},
		{"threads", required_argument, 0, 0},
		{"kmer", required_argument, 0, 0},
		{"min-ovlp", required_argument, 0, 0},
		{"extra-params", required_argument, 0, 0},
		{"meta", no_argument, 0, 0},
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
			else if (!strcmp(longOptions[optionIndex].name, "min-read"))
				minReadLength = atoi(optarg);
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
			else if (!strcmp(longOptions[optionIndex].name, "reads"))
				readsFasta = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "out-asm"))
				outAssembly = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "genome-size"))
				genomeSize = atoll(optarg);
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
	if (readsFasta.empty() || outAssembly.empty() || 
		configPath.empty())
	{
		printUsage();
		return false;
	}

	std::cout << "readsFasta: " << readsFasta << std::endl;
	std::cout << "outAssembly: " << outAssembly << std::endl;
	std::cout << "configPath: " << configPath << std::endl;

	return true;
}

int assemble_main(int argc, char** argv)
{
	#ifdef NDEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	int kmerSize = -1;
	int minReadLength = 0;
	size_t genomeSize = 0;
	int minOverlap = 5000;
	bool debugging = false;
	bool unevenCov = false;
	size_t numThreads = 1;
	std::string readsFasta;
	std::string outAssembly;
	std::string logFile;
	std::string configPath;
	std::string extraParams;

	if (!parseArgs(argc, argv, readsFasta, outAssembly, logFile, genomeSize,
				   kmerSize, debugging, numThreads, minOverlap, configPath, 
				   minReadLength, unevenCov, extraParams)) return 1;

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
	Parameters::get().unevenCoverage = unevenCov;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 
	Logger::get().debug() << "Running with minimum overlap " << minOverlap;
	Logger::get().debug() << "Metagenome mode: " << "NY"[unevenCov];

	//TODO: unify minimumOverlap ad safeOverlap concepts
	Parameters::get().minimumOverlap = 1000;

	SequenceContainer readsContainer;
	std::vector<std::string> readsList;
	readsList.push_back(readsFasta);
	Logger::get().info() << "Reading sequences";
	try
	{
		//only use reads that are longer than minOverlap,
		//or a specified threshold (used for downsampling)
		minReadLength = std::max(minReadLength, minOverlap);
		for (auto& readsFile : readsList)
		{
			readsContainer.loadFromFile(readsFile, minReadLength);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	readsContainer.buildPositionIndex();
	VertexIndex vertexIndex(readsContainer, 
							(int)Config::get("assemble_kmer_sample"));
	vertexIndex.outputProgress(true);

	/*int64_t sumLength = 0;
	for (auto& seq : readsContainer.iterSeqs())
	{
		sumLength += seq.sequence.length();
	}
	int coverage = sumLength / 2 / genomeSize;
	Logger::get().debug() << "Expected read coverage: " << coverage;*/

	const int MIN_FREQ = 2;
	static const float SELECT_RATE = Config::get("meta_read_top_kmer_rate");
	static const int TANDEM_FREQ = Config::get("meta_read_filter_kmer_freq");

	//Building index
	bool useMinimizers = Config::get("use_minimizers");
	if (useMinimizers)
	{
		const int minWnd = Config::get("minimizer_window");
		vertexIndex.buildIndexMinimizers(/*min freq*/ 1, minWnd);
	}
	else	//indexing using solid k-mers
	{
		vertexIndex.countKmers();
		vertexIndex.buildIndexUnevenCoverage(MIN_FREQ, SELECT_RATE, 
											 TANDEM_FREQ);
	}

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	//int maxOverlapsNum = !Parameters::get().unevenCoverage ? 5 * coverage : 0;
	OverlapDetector ovlp(readsContainer, vertexIndex,
						 (int)Config::get("maximum_jump"), 
						 Parameters::get().minimumOverlap,
						 (int)Config::get("maximum_overhang"),
						 /*store alignment*/ false,
						 /*only max ovlp*/ true,
						 /*no div threshold*/ 1.0f,
						 (bool)Config::get("reads_base_alignment"),
						 /*partition bad map*/ false,
						 (bool)Config::get("hpc_scoring_on"));
	OverlapContainer readOverlaps(ovlp, readsContainer);
	readOverlaps.estimateOverlaperParameters();
	readOverlaps.setDivergenceThreshold((float)Config::get("assemble_ovlp_divergence"),
										(bool)Config::get("assemble_divergence_relative"));

	Extender extender(readsContainer, readOverlaps, minOverlap);
	extender.assembleDisjointigs();
	vertexIndex.clear();

	ConsensusGenerator consGen;
	auto disjointigsFasta = 
		consGen.generateConsensuses(extender.getDisjointigPaths());
	SequenceContainer::writeFasta(disjointigsFasta, outAssembly);

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	return 0;
}
