#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Quick and dirty alignment consensus
"""

from __future__ import absolute_import
from __future__ import division
import logging
from collections import defaultdict
from flye.six.moves import range
from flye.six import itervalues

import multiprocessing
import signal
import traceback

from flye.polishing.alignment import shift_gaps, get_uniform_alignments
from flye.utils.sam_parser import SynchronizedSamReader
import flye.config.py_cfg as cfg
import flye.utils.fasta_parser as fp
from flye.six.moves import zip

logger = logging.getLogger()

class Profile(object):
    __slots__ = ("insertions", "matches", "nucl")

    def __init__(self):
        self.insertions = defaultdict(str)
        self.matches = defaultdict(int)
        self.nucl = "-"

def _thread_worker(aln_reader, contigs_info, platform, results_queue,
                   error_queue):
    try:
        while not aln_reader.is_eof():
            ctg_id, ctg_aln = aln_reader.get_chunk()
            if ctg_id is None:
                break

            profile, aln_errors = _contig_profile(ctg_aln, platform,
                                                  contigs_info[ctg_id].length)
            sequence = _flatten_profile(profile)
            results_queue.put((ctg_id, sequence, aln_errors))

    except Exception as e:
        logger.error("Thread exception")
        logger.error(traceback.format_exc())
        error_queue.put(e)


def get_consensus(alignment_path, contigs_path, contigs_info, num_proc,
                  platform):
    """
    Main function
    """
    logger.info("get_consensus(alignment_path: {}, contigs_path: {})".format(alignment_path, contigs_path))
    aln_reader = SynchronizedSamReader(alignment_path,
                                       fp.read_sequence_dict(contigs_path),
                                       max_coverage=cfg.vals["max_read_coverage"],
                                       use_secondary=True)
    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    error_queue = manager.Queue()

    #making sure the main process catches SIGINT
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    threads = []
    for _ in range(num_proc):
        threads.append(multiprocessing.Process(target=_thread_worker,
                                               args=(aln_reader, contigs_info,
                                                     platform, results_queue,
                                                     error_queue)))
    signal.signal(signal.SIGINT, orig_sigint)

    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
            if t.exitcode == -9:
                logger.error("Looks like the system ran out of memory")
            if t.exitcode != 0:
                raise Exception("One of the processes exited with code: {0}"
                                .format(t.exitcode))
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()
        raise

    if not error_queue.empty():
        raise error_queue.get()
    aln_reader.close()

    out_fasta = {}
    total_aln_errors = []
    while not results_queue.empty():
        ctg_id, ctg_seq, aln_errors = results_queue.get()
        total_aln_errors.extend(aln_errors)
        if len(ctg_seq) > 0:
            out_fasta[ctg_id] = ctg_seq

    mean_aln_error = sum(total_aln_errors) / (len(total_aln_errors) + 1)
    logger.info("Alignment error rate: %f", mean_aln_error)

    return out_fasta


def _contig_profile(alignment, platform, genome_len):
    """
    Computes alignment profile
    """

    #leave the best uniform alignments
    alignment = get_uniform_alignments(alignment, genome_len)

    aln_errors = []
    profile = [Profile() for _ in range(genome_len)]
    #max_aln_err = cfg.vals["err_modes"][platform]["max_aln_error"]
    for aln in alignment:
        #if aln.err_rate > max_aln_err: continue
        aln_errors.append(aln.err_rate)

        #after gap shifting it is possible that
        #two gaps are aligned against each other
        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in zip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            if trg_pos >= genome_len:
                trg_pos -= genome_len

            #total += 1
            prof_elem = profile[trg_pos]
            if trg_nuc == "-" and qry_nuc != "-":
                prof_elem.insertions[aln.qry_id] += qry_nuc
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.matches[qry_nuc] += 1

            trg_pos += 1

    #print "len", genome_len, "median coverage", cov_threshold
    #print "total bases: ", total, "discarded bases: ", discarded
    #print "filtered", float(discarded) / total
    #print ""

    return profile, aln_errors


def _flatten_profile(profile):
    growing_seq = []
    ins_group = defaultdict(int)

    for elem in profile:
        pos_matches = elem.matches
        pos_insertions = elem.insertions
        pos_nucl = elem.nucl

        ins_group.clear()
        for ins_str in itervalues(pos_insertions):
            ins_group[ins_str] += 1

        match_and_del_num = sum(itervalues(pos_matches))
        del_num = pos_matches["-"]
        num_ins = len(pos_insertions)

        max_match = pos_nucl
        if len(pos_matches):
            max_match = max(sorted(pos_matches), key=pos_matches.get)
        max_insert = None
        if ins_group:
            max_insert = max(sorted(ins_group), key=ins_group.get)

        is_deletion = max_match == "-" or del_num > match_and_del_num // 3
        is_insertion = max_insert and max_insert != "-" and num_ins > match_and_del_num // 3

        if not is_deletion:
            growing_seq.append(max_match)
        if is_insertion:
            growing_seq.append(max_insert)

    return "".join(growing_seq)
