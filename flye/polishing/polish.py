#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs polishing binary in parallel and concatentes output
"""

import logging
import random
import subprocess
import os
from collections import defaultdict
from threading import Thread

from flye.polishing.alignment import (make_alignment, get_contigs_info,
                                      SynchronizedSamReader, merge_chunks,
                                      split_into_chunks)
from flye.polishing.bubbles import make_bubbles
import flye.utils.fasta_parser as fp
from flye.utils.utils import which
import flye.config.py_cfg as cfg


POLISH_BIN = "flye-polish"

logger = logging.getLogger()


class PolishException(Exception):
    pass


def check_binaries():
    if not which(POLISH_BIN):
        raise PolishException("polishing binary was not found. "
                              "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([POLISH_BIN, "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))
    except OSError as e:
        raise PolishException(str(e))


def polish(contig_seqs, read_seqs, work_dir, num_iters, num_threads, error_mode,
           output_progress):
    """
    High-level polisher interface
    """
    logger_state = logger.disabled
    if not output_progress:
        logger.disabled = True

    subs_matrix = os.path.join(cfg.vals["pkg_root"],
                               cfg.vals["err_modes"][error_mode]["subs_matrix"])
    hopo_matrix = os.path.join(cfg.vals["pkg_root"],
                               cfg.vals["err_modes"][error_mode]["hopo_matrix"])
    stats_file = os.path.join(work_dir, "contigs_stats.txt")

    prev_assembly = contig_seqs
    contig_lengths = None
    coverage_stats = None
    for i in xrange(num_iters):
        logger.info("Polishing genome ({0}/{1})".format(i + 1, num_iters))

        #split into 1Mb chunks to reduce RAM usage
        #slightly vary chunk size between iterations
        CHUNK_SIZE = 1000000 - (i % 2) * 100000
        chunks_file = os.path.join(work_dir, "chunks_{0}.fasta".format(i + 1))
        chunks = split_into_chunks(fp.read_sequence_dict(prev_assembly),
                                       CHUNK_SIZE)
        fp.write_fasta_dict(chunks, chunks_file)

        ####
        logger.info("Running minimap2")
        alignment_file = os.path.join(work_dir, "minimap_{0}.sam".format(i + 1))
        make_alignment(chunks_file, read_seqs, num_threads,
                       work_dir, error_mode, alignment_file,
                       reference_mode=True, sam_output=True)

        #####
        logger.info("Separating alignment into bubbles")
        contigs_info = get_contigs_info(chunks_file)
        bubbles_file = os.path.join(work_dir,
                                    "bubbles_{0}.fasta".format(i + 1))
        coverage_stats, mean_aln_error = \
            make_bubbles(alignment_file, contigs_info, chunks_file,
                         error_mode, num_threads,
                         bubbles_file)

        logger.info("Alignment error rate: {0}".format(mean_aln_error))
        consensus_out = os.path.join(work_dir, "consensus_{0}.fasta".format(i + 1))
        polished_file = os.path.join(work_dir, "polished_{0}.fasta".format(i + 1))
        if os.path.getsize(bubbles_file) == 0:
            logger.info("No reads were aligned during polishing")
            if not output_progress:
                logger.disabled = logger_state
            open(stats_file, "w").write("seq_name\tlength\tcoverage\n")
            open(polished_file, "w")
            return polished_file, stats_file

        #####
        logger.info("Correcting bubbles")
        _run_polish_bin(bubbles_file, subs_matrix, hopo_matrix,
                        consensus_out, num_threads, output_progress)
        polished_fasta, polished_lengths = _compose_sequence(consensus_out)
        merged_chunks = merge_chunks(polished_fasta)
        fp.write_fasta_dict(merged_chunks, polished_file)

        #Cleanup
        os.remove(chunks_file)
        os.remove(bubbles_file)
        os.remove(consensus_out)
        os.remove(alignment_file)

        contig_lengths = polished_lengths
        prev_assembly = polished_file

    #merge information from chunks
    contig_lengths = merge_chunks(contig_lengths, fold_function=sum)
    coverage_stats = merge_chunks(coverage_stats,
                                  fold_function=lambda l: sum(l) / len(l))

    with open(stats_file, "w") as f:
        f.write("seq_name\tlength\tcoverage\n")
        for ctg_id in contig_lengths:
            f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                    contig_lengths[ctg_id], coverage_stats[ctg_id]))

    if not output_progress:
        logger.disabled = logger_state

    return prev_assembly, stats_file


def generate_polished_edges(edges_file, gfa_file, polished_contigs, work_dir,
                            error_mode, num_threads):
    """
    Generate polished graph edges sequences by extracting them from
    polished contigs
    """
    logger.debug("Generating polished GFA")

    alignment_file = os.path.join(work_dir, "edges_aln.sam")
    polished_dict = fp.read_sequence_dict(polished_contigs)
    make_alignment(polished_contigs, [edges_file], num_threads,
                   work_dir, error_mode, alignment_file,
                   reference_mode=True, sam_output=True)
    aln_reader = SynchronizedSamReader(alignment_file,
                                       polished_dict,
                                       cfg.vals["max_read_coverage"])
    aln_reader.init_reading()
    aln_by_edge = defaultdict(list)

    #getting one best alignment for each contig
    while not aln_reader.is_eof():
        _, ctg_aln = aln_reader.get_chunk()
        for aln in ctg_aln:
            aln_by_edge[aln.qry_id].append(aln)

    MIN_CONTAINMENT = 0.9
    updated_seqs = 0
    edges_dict = fp.read_sequence_dict(edges_file)
    for edge in edges_dict:
        if edge in aln_by_edge:
            main_aln = aln_by_edge[edge][0]
            map_start = main_aln.trg_start
            map_end = main_aln.trg_end
            for aln in aln_by_edge[edge]:
                if aln.trg_id == main_aln.trg_id and aln.trg_sign == main_aln.trg_sign:
                    map_start = min(map_start, aln.trg_start)
                    map_end = max(map_end, aln.trg_end)

            new_seq = polished_dict[main_aln.trg_id][map_start : map_end]
            if main_aln.qry_sign == "-":
                new_seq = fp.reverse_complement(new_seq)

            #print edge, main_aln.qry_len, len(new_seq), main_aln.qry_start, main_aln.qry_end
            if float(len(new_seq)) / aln.qry_len > MIN_CONTAINMENT:
                edges_dict[edge] = new_seq
                updated_seqs += 1

    #writes fasta file with polished egdes
    edges_polished = os.path.join(work_dir, "polished_edges.fasta")
    fp.write_fasta_dict(edges_dict, edges_polished)

    #writes gfa file with polished edges
    gfa_polished = open(os.path.join(work_dir, "polished_edges.gfa"), "w")
    for line in open(gfa_file, "r"):
        if line.startswith("S"):
            seq_id = line.split()[1]
            coverage_tag = line.split()[3]
            gfa_polished.write("S\t{0}\t{1}\t{2}\n"
                                .format(seq_id, edges_dict[seq_id], coverage_tag))
        else:
            gfa_polished.write(line)

    logger.debug("{0} sequences remained unpolished"
                    .format(len(edges_dict) - updated_seqs))
    os.remove(alignment_file)


def _run_polish_bin(bubbles_in, subs_matrix, hopo_matrix,
                    consensus_out, num_threads, output_progress):
    """
    Invokes polishing binary
    """
    cmdline = [POLISH_BIN, "-t", str(num_threads)]
    if not output_progress:
        cmdline.append("-q")
    cmdline.extend([bubbles_in, subs_matrix,
                    hopo_matrix, consensus_out])

    try:
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))
    except OSError as e:
        raise PolishException(str(e))


def _compose_sequence(consensus_file):
    """
    Concatenates bubbles consensuses into genome
    """
    consensuses = defaultdict(list)
    coverage = defaultdict(list)
    with open(consensus_file, "r") as f:
        header = True
        for line in f:
            if header:
                tokens = line.strip().split(" ")
                ctg_id = tokens[0][1:]
                ctg_pos = int(tokens[1])
                coverage[ctg_id].append(int(tokens[2]))
            else:
                consensuses[ctg_id].append((ctg_pos, line.strip()))
            header = not header

    polished_fasta = {}
    polished_stats = {}
    for ctg_id, seqs in consensuses.iteritems():
        sorted_seqs = map(lambda p: p[1], sorted(seqs, key=lambda p: p[0]))
        concat_seq = "".join(sorted_seqs)
        mean_coverage = sum(coverage[ctg_id]) / len(coverage[ctg_id])
        polished_fasta[ctg_id] = concat_seq
        polished_stats[ctg_id] = len(concat_seq)

    return polished_fasta, polished_stats
