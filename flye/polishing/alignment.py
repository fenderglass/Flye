#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs Minimap2 and parses its output
"""

from __future__ import absolute_import
from __future__ import division
import os
from collections import namedtuple
import subprocess
import logging
import datetime

import flye.utils.fasta_parser as fp
from flye.utils.utils import which
from flye.utils.sam_parser import AlignmentException
from flye.six import iteritems
from flye.six.moves import range


logger = logging.getLogger()
MINIMAP_BIN = "flye-minimap2"
SAMTOOLS_BIN = "flye-samtools"

ContigInfo = namedtuple("ContigInfo", ["id", "length", "type"])


def check_binaries():
    if not which(MINIMAP_BIN):
        raise AlignmentException("minimap2 is not installed")
    if not which(SAMTOOLS_BIN):
        raise AlignmentException("samtools is not installed")
    if not which("sort"):
        raise AlignmentException("UNIX sort utility is not available")


def make_alignment(reference_file, reads_file, num_proc,
                   work_dir, platform, out_alignment, reference_mode,
                   sam_output):
    """
    Runs minimap2 and sorts its output
    """
    minimap_ref_mode = {False: "ava", True: "map"}
    minimap_reads_mode = {"nano": "ont", "pacbio": "pb"}
    mode = minimap_ref_mode[reference_mode] + "-" + minimap_reads_mode[platform]

    _run_minimap(reference_file, reads_file, num_proc, mode,
                 out_alignment, sam_output)

    #if sam_output:
    #    preprocess_sam(out_alignment, work_dir)


def get_contigs_info(contigs_file):
    contigs_info = {}
    contigs_fasta = fp.read_sequence_dict(contigs_file)
    for ctg_id, ctg_seq in iteritems(contigs_fasta):
        contig_type = ctg_id.split("_")[0]
        contigs_info[ctg_id] = ContigInfo(ctg_id, len(ctg_seq),
                                          contig_type)

    return contigs_info


def shift_gaps(seq_trg, seq_qry):
    """
    Shifts all ambigious query gaps to the right
    """
    lst_trg, lst_qry = list("$" + seq_trg + "$"), list("$" + seq_qry + "$")
    is_gap = False
    gap_start = 0
    for i in range(len(lst_trg)):
        if is_gap and lst_qry[i] != "-":
            is_gap = False
            swap_left = gap_start - 1
            swap_right = i - 1

            while (swap_left > 0 and swap_right >= gap_start and
                   lst_qry[swap_left] == lst_trg[swap_right]):
                lst_qry[swap_left], lst_qry[swap_right] = \
                            lst_qry[swap_right], lst_qry[swap_left]
                swap_left -= 1
                swap_right -= 1

        if not is_gap and lst_qry[i] == "-":
            is_gap = True
            gap_start = i

    return "".join(lst_qry[1 : -1])


def get_uniform_alignments(alignments, seq_len):
    """
    Leaves top alignments for each position within contig
    assuming uniform coverage distribution
    """
    def _get_median(lst):
        if not lst:
            raise ValueError("_get_median() arg is an empty sequence")
        sorted_list = sorted(lst)
        if len(lst) % 2 == 1:
            return sorted_list[len(lst) // 2]
        else:
            mid1 = sorted_list[(len(lst) // 2) - 1]
            mid2 = sorted_list[(len(lst) // 2)]
            return (mid1 + mid2) / 2

    WINDOW = 100
    MIN_COV = 10
    COV_RATE = 1.25

    #split contig into windows, get median read coverage over all windows and
    #determine the quality threshold cutoffs for each window
    wnd_primary_cov = [0 for _ in range(seq_len // WINDOW + 1)]
    wnd_aln_quality = [[] for _ in range(seq_len // WINDOW + 1)]
    wnd_qual_thresholds = [1.0 for _ in range(seq_len // WINDOW + 1)]
    for aln in alignments:
        for i in range(aln.trg_start // WINDOW, aln.trg_end // WINDOW):
            if not aln.is_secondary:
                wnd_primary_cov[i] += 1
            wnd_aln_quality[i].append(aln.err_rate)

    #for each window, select top X alignmetns, where X is the median read coverage
    cov_threshold = max(int(COV_RATE * _get_median(wnd_primary_cov)), MIN_COV)
    for i in range(len(wnd_aln_quality)):
        if len(wnd_aln_quality[i]) > cov_threshold:
            wnd_qual_thresholds[i] = sorted(wnd_aln_quality[i])[cov_threshold]

    #for each alignment, count in how many windows it passes the threshold
    filtered_alignments = []
    total_sequence = 0
    filtered_sequence = 0
    for aln in alignments:
        good_windows = 0
        total_windows = aln.trg_end // WINDOW - aln.trg_start // WINDOW
        total_sequence += aln.trg_end - aln.trg_start
        for i in range(aln.trg_start // WINDOW, aln.trg_end // WINDOW):
            if aln.err_rate <= wnd_qual_thresholds[i]:
                good_windows += 1

        if good_windows > total_windows // 2:
            filtered_alignments.append(aln)
            filtered_sequence += aln.trg_end - aln.trg_start

    #filtered_reads_rate = 1 - float(len(filtered_alignments)) / len(alignments)
    #filtered_seq_rate = 1 - float(filtered_sequence) / total_sequence
    #logger.debug("Filtered {0:7.2f}% reads, {1:7.2f}% sequence"
    #                .format(filtered_reads_rate * 100, filtered_seq_rate * 100))

    return filtered_alignments


def split_into_chunks(fasta_in, chunk_size):
    out_dict = {}
    for header, seq in iteritems(fasta_in):
        #print len(seq)
        for i in range(0, max(len(seq) // chunk_size, 1)):
            chunk_hdr = "{0}$chunk_{1}".format(header, i)
            start = i * chunk_size
            end = (i + 1) * chunk_size
            if len(seq) - end < chunk_size:
                end = len(seq)

            #print(start, end)
            out_dict[chunk_hdr] = seq[start : end]

    return out_dict


def merge_chunks(fasta_in, fold_function=lambda l: "".join(l)):
    """
    Merges sequence chunks. Chunk names are in format `orig_name$chunk_id`.
    Each chunk is as dictionary entry. Value type is arbitrary and
    one can supply a custom fold function
    """
    def name_split(h):
        orig_hdr, chunk_id = h.rsplit("$", 1)
        return orig_hdr, int(chunk_id.rsplit("_", 1)[1])

    out_dict = {}
    cur_seq = []
    cur_contig = None
    for hdr in sorted(fasta_in, key=name_split):
        orig_name, dummy_chunk_id = name_split(hdr)
        if orig_name != cur_contig:
            if cur_contig != None:
                out_dict[cur_contig] = fold_function(cur_seq)
            cur_seq = []
            cur_contig = orig_name
        cur_seq.append(fasta_in[hdr])

    if cur_seq:
        out_dict[cur_contig] = fold_function(cur_seq)

    return out_dict


def _run_minimap(reference_file, reads_files, num_proc, mode, out_file,
                 sam_output):
    #SAM_HEADER = "\'@PG|@HD|@SQ|@RG|@CO\'"
    work_dir = os.path.dirname(out_file)
    stderr_file = os.path.join(work_dir, "minimap.stderr")
    SORT_THREADS = "4"
    SORT_MEM = "4G" if os.path.getsize(reference_file) > 100 * 1024 * 1024 else "1G"

    cmdline = [MINIMAP_BIN, "'" + reference_file + "'"]
    cmdline.extend(["'" + read_file + "'" for read_file in reads_files])
    cmdline.extend(["-x", mode, "-t", str(num_proc)])

    #Produces gzipped SAM sorted by reference name. Since it's not sorted by
    #read name anymore, it's important that all reads have SEQ.
    #is sam_output not set, produces PAF alignment
    #a = SAM output, p = min primary-to-seconday score
    #N = max secondary alignments
    #--sam-hit-only = don't output unmapped reads
    #--secondary-seq = custom option to output SEQ for seqcondary alignment with hard clipping
    #-L: move CIGAR strings for ultra-long reads to the separate tag
    #-Q don't output fastq quality
    if sam_output:
        tmp_prefix = os.path.join(os.path.dirname(out_file),
                                  "sort_" + datetime.datetime.now().strftime("%y%m%d_%H%M%S"))
        cmdline.extend(["-a", "-p", "0.5", "-N", "10", "--sam-hit-only", "-L",
                        "-z", "1000", "-Q", "--secondary-seq", "-I", "64G"])
        cmdline.extend(["|", SAMTOOLS_BIN, "view", "-T", "'" + reference_file + "'", "-u", "-"])
        cmdline.extend(["|", SAMTOOLS_BIN, "sort", "-T", "'" + tmp_prefix + "'", "-O", "bam",
                        "-@", SORT_THREADS, "-l", "1", "-m", SORT_MEM])
    else:
        pass    #paf output enabled by default

        #cmdline.extend(["|", "grep", "-Ev", SAM_HEADER])    #removes headers
        #cmdline.extend(["|", "sort", "-k", "3,3", "-T", work_dir,
        #                "--parallel=8", "-S", "4G"])
        #cmdline.extend(["|", "gzip", "-1"])

    #logger.debug("Running: " + " ".join(cmdline))
    try:
        devnull = open(os.devnull, "wb")
        #env = os.environ.copy()
        #env["LC_ALL"] = "C"

        command_line = ["/bin/bash", "-c",
                              "set -eo pipefail; " + " ".join(cmdline)]
        logger.info("Running minimap2 with the following command line: {}".format(command_line))
        subprocess.check_call(command_line,
                              stderr=open(stderr_file, "w"),
                              stdout=open(out_file, "w"))
        os.remove(stderr_file)

    except (subprocess.CalledProcessError, OSError) as e:
        logger.error("Error running minimap2, terminating. See the alignment error log "
                     " for details: " + stderr_file)
        raise AlignmentException(str(e))
