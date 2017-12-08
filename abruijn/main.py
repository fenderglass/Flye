#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Main logic of the package
"""

from __future__ import print_function
import sys
import os
import logging
import argparse
import json

import abruijn.alignment as aln
import abruijn.bubbles as bbl
import abruijn.polish as pol
import abruijn.fasta_parser as fp
import abruijn.assemble as asm
import abruijn.repeat_graph as repeat
import abruijn.consensus as cons
from abruijn.__version__ import __version__
import abruijn.config as config


logger = logging.getLogger()

class ResumeException(Exception):
    pass

class Job(object):
    """
    Describes an abstract list of jobs with persistent
    status that can be resumed
    """
    run_description = {"stage_name" : ""}

    def __init__(self):
        self.name = None
        self.args = None
        self.work_dir = None
        self.out_files = []

    def run(self):
        pass

    def save(self, save_file):
        Job.run_description["stage_name"] = self.name

        with open(save_file, "w") as fp:
            json.dump(Job.run_description, fp)

    def load(self, save_file):
        with open(save_file, "r") as fp:
            data = json.load(fp)
            Job.run_description = data

    def completed(self, save_file):
        with open(save_file, "r") as fp:
            data = json.load(fp)

            for file in self.out_files:
                if not os.path.exists(file):
                    return False

            return True


class JobAssembly(Job):
    def __init__(self, out_assembly, log_file):
        super(JobAssembly, self).__init__()
        self.out_assembly = out_assembly
        self.log_file = log_file
        self.name = "assembly"
        self.out_files = [out_assembly]

    def run(self):
        asm.assemble(self.args, self.out_assembly, self.log_file)


class JobRepeat(Job):
    def __init__(self, in_assembly, out_folder, log_file):
        super(JobRepeat, self).__init__()
        self.in_assembly = in_assembly
        self.log_file = log_file
        self.out_folder = out_folder
        self.name = "repeat"

        #edges_sequences = os.path.join(out_folder, "graph_final.fasta")
        edges_sequences = os.path.join(out_folder, "graph_paths.fasta")
        repeat_graph = os.path.join(out_folder, "graph_final.gfa")
        self.out_files = [edges_sequences, repeat_graph]

    def run(self):
        logger.info("Performing repeat analysis")
        repeat.analyse_repeats(self.args, self.in_assembly, self.out_folder,
                               self.log_file)


class JobAlignment(Job):
    def __init__(self, in_reference, out_alignment, stage_id):
        super(JobAlignment, self).__init__()
        self.in_reference = in_reference
        self.out_alignment = out_alignment
        self.name = "alignment_" + str(stage_id)
        self.stage_id = stage_id
        self.out_files = [out_alignment]

    def run(self):
        #logger.info("Polishing genome ({0}/{1})".format(self.stage_id,
        #                                                self.args.num_iters))
        logger.info("Running Minimap2")
        aln.make_alignment(self.in_reference, self.args.reads, self.args.threads,
                           self.work_dir, self.args.platform, self.out_alignment)


class JobConsensus(Job):
    def __init__(self, in_contigs, in_alignment, out_consensus):
        super(JobConsensus, self).__init__()

        self.in_contigs = in_contigs
        self.in_alignment = in_alignment
        self.out_consensus = out_consensus
        self.name = "consensus"
        self.out_files = [out_consensus]

    def run(self):
        logger.info("Computing rough consensus")
        contigs_info = aln.get_contigs_info(self.in_contigs)
        consensus_fasta = cons.get_consensus(self.in_alignment, self.in_contigs,
                                             contigs_info, self.args.threads,
                                             self.args.platform,
                                             self.args.min_overlap)
        fp.write_fasta_dict(consensus_fasta, self.out_consensus)


class JobPolishing(Job):
    def __init__(self, in_contigs, in_alignment, out_consensus,
                 stage_id, seq_platform):
        super(JobPolishing, self).__init__()

        self.in_contigs = in_contigs
        self.in_alignment = in_alignment
        self.out_consensus = out_consensus
        self.name = "polishing_" + str(stage_id)
        self.stage_id = stage_id
        self.seq_platform = seq_platform
        self.out_files = [out_consensus]

    def run(self):
        logger.info("Polishing genome ({0}/{1})".format(self.stage_id,
                                                        self.args.num_iters))
        contigs_info = aln.get_contigs_info(self.in_contigs)

        logger.info("Separating alignment into bubbles")
        bubbles = bbl.get_bubbles(self.in_alignment, contigs_info, self.in_contigs,
                                  self.seq_platform, self.args.threads,
                                  self.args.min_overlap)
        logger.info("Correcting bubbles")
        polished_fasta = pol.polish(bubbles, self.args.threads,
                                    self.seq_platform, self.work_dir,
                                    self.stage_id)
        fp.write_fasta_dict(polished_fasta, self.out_consensus)


def _create_job_list(args, work_dir, log_file):
    """
    Build pipeline as a list of consecutive jobs
    """
    jobs = []

    #Assembly job
    draft_assembly = os.path.join(work_dir, "draft_assembly.fasta")
    jobs.append(JobAssembly(draft_assembly, log_file))

    #Pre-polishing
    alignment_file = os.path.join(work_dir, "minimap_0.sam")
    pre_polished_file = os.path.join(work_dir, "polished_0.fasta")
    jobs.append(JobAlignment(draft_assembly, alignment_file, 0))
    jobs.append(JobConsensus(draft_assembly, alignment_file,
                             pre_polished_file))

    #Repeat analysis
    #edges_sequences = os.path.join(work_dir, "graph_final.fasta")
    edges_sequences = os.path.join(work_dir, "graph_paths.fasta")
    jobs.append(JobRepeat(pre_polished_file, work_dir, log_file))

    #Full polishing
    prev_assembly = edges_sequences
    for i in xrange(args.num_iters):
        alignment_file = os.path.join(work_dir, "minimap_{0}.sam".format(i + 1))
        polished_file = os.path.join(work_dir,
                                     "polished_{0}.fasta".format(i + 1))
        jobs.append(JobAlignment(prev_assembly, alignment_file, i + 1))
        jobs.append(JobPolishing(prev_assembly, alignment_file, polished_file,
                                 i + 1, args.platform))
        prev_assembly = polished_file

    for job in jobs:
        job.args = args
        job.work_dir = work_dir

    return jobs


def _get_kmer_size(args):
    """
    Select k-mer size based on the target genome size
    """
    suffix = args.reads.rsplit(".", 1)[-1]
    if suffix in ["fasta", "fa"]:
        reads_size = os.path.getsize(args.reads)
    elif suffix in ["fastq", "fq"]:
        reads_size = os.path.getsize(args.reads) / 2
    else:
        raise ResumeException("Uknown input reads format: " + suffix)

    genome_size = reads_size / args.coverage
    logger.debug("Estimated genome size: {0}".format(genome_size))
    kmer_size = 15
    if genome_size > config.vals["big_genome"]:
        kmer_size = 17
    logger.debug("Chosen k-mer size: {0}".format(kmer_size))
    return kmer_size


def _run(args):
    """
    Runs the pipeline
    """
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    work_dir = os.path.abspath(args.out_dir)

    log_file = os.path.join(work_dir, "abruijn.log")
    _enable_logging(log_file, args.debug,
                    overwrite=not args.resume and not args.resume_from)

    logger.info("Running ABruijn")
    logger.info("  Run invoked with following parameters:")
    for arg in vars(args):
        logger.info("    {}: {}".format(arg,getattr(args, arg)))    
    aln.check_binaries()
    pol.check_binaries()
    asm.check_binaries()

    if not os.path.exists(args.reads):
        raise ResumeException("Can't open " + args.reads)

    if args.kmer_size is None:
        args.kmer_size = _get_kmer_size(args)

    save_file = os.path.join(work_dir, "abruijn.save")
    jobs = _create_job_list(args, work_dir, log_file)

    current_job = 0
    if args.resume or args.resume_from:
        if not os.path.exists(save_file):
            raise ResumeException("Can't find save file")

        logger.info("Resuming previous run")
        if args.resume_from:
            job_to_resume = args.resume_from
        else:
            job_to_resume = json.load(open(save_file, "r"))["stage_name"]

        for i in xrange(len(jobs)):
            if jobs[i].name == job_to_resume:
                jobs[i].load(save_file)
                current_job = i
                if not jobs[i - 1].completed(save_file):
                    raise ResumeException("Can't resume: stage {0} incomplete"
                                          .format(jobs[i].name))
                break

    for i in xrange(current_job, len(jobs)):
        jobs[i].save(save_file)
        jobs[i].run()

    logger.info("Done! Your assembly is in file: " + jobs[-1].out_files[0])


def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def main():
    def check_int_range(value, min_val, max_val, require_odd=False):
        ival = int(value)
        if ival < min_val or ival > max_val:
             raise argparse.ArgumentTypeError("value should be in "
                            "range [{0}, {1}]".format(min_val, max_val))
        if require_odd and ival % 2 == 0:
            raise argparse.ArgumentTypeError("should be an odd number")
        return ival

    parser = argparse.ArgumentParser(description="ABruijn: assembly of long and"
                                     " error-prone reads")

    parser.add_argument("reads", metavar="reads",
                        help="path to reads file (FASTA/Q format)")
    parser.add_argument("out_dir", metavar="out_dir",
                        help="output directory")
    parser.add_argument("coverage", metavar="coverage (integer)",
                        type=lambda v: check_int_range(v, 1, 1000),
                        help="estimated assembly coverage")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("--resume", action="store_true",
                        dest="resume", default=False,
                        help="resume from the last completed stage")
    parser.add_argument("--resume-from", dest="resume_from",
                        default=None, help="resume from a custom stage")
    parser.add_argument("-t", "--threads", dest="threads",
                        type=lambda v: check_int_range(v, 1, 128),
                        default=1, help="number of parallel threads "
                        "(default: 1)")
    parser.add_argument("-i", "--iterations", dest="num_iters",
                        type=lambda v: check_int_range(v, 0, 10),
                        default=1, help="number of polishing iterations "
                        "(default: 1)")
    parser.add_argument("-p", "--platform", dest="platform",
                        default="pacbio",
                        choices=["pacbio", "nano", "pacbio_hi_err"],
                        help="sequencing platform (default: pacbio)")
    parser.add_argument("-k", "--kmer-size", dest="kmer_size",
                        type=lambda v: check_int_range(v, 11, 31, require_odd=True),
                        default=None, help="kmer size (default: auto)")
    parser.add_argument("-o", "--min-overlap", dest="min_overlap",
                        type=lambda v: check_int_range(v, 2000, 10000),
                        default=5000, help="minimum overlap between reads "
                        "(default: 5000)")
    parser.add_argument("-m", "--min-coverage", dest="min_kmer_count",
                        type=lambda v: check_int_range(v, 1, 1000),
                        default=None, help="minimum kmer coverage "
                        "(default: auto)")
    parser.add_argument("-x", "--max-coverage", dest="max_kmer_count",
                        type=lambda v: check_int_range(v, 1, 1000),
                        default=None, help="maximum kmer coverage "
                        "(default: auto)")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    try:
        _run(args)
    except (aln.AlignmentException, pol.PolishException,
            asm.AssembleException, repeat.RepeatException,
            ResumeException) as e:
        logger.error(e)
        return 1

    return 0
