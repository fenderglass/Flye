Flye manual
===========

Table of Contents
-----------------

- [Quick usage](#quickusage)
- [Examples](#examples)
- [Supported Input Data](#inputdata)
- [Parameter Descriptions](#parameters)
- [Flye output](#output)
- [Repeat graph](#graph)
- [Flye benchmarks](#performance)
- [Algorithm Description](#algorithm)


## <a name="quickusage"></a> Quick usage

```
usage: flye (--pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw |
	     --nano-corr | --subassemblies) file1 [file_2 ...]
	     --out-dir PATH

	     [--genome-size SIZE] [--threads int] [--iterations int]
	     [--meta] [--plasmids] [--trestle] [--polish-target]
	     [--keep-haplotypes] [--debug] [--version] [--help] 
	     [--resume] [--resume-from] [--stop-after] [--min-overlap SIZE]

Assembly of long reads with repeat graphs

optional arguments:
  -h, --help            show this help message and exit
  --pacbio-raw path [path ...]
                        PacBio raw reads
  --pacbio-corr path [path ...]
                        PacBio corrected reads
  --pacbio-hifi path [path ...]
                        PacBio HiFi reads
  --nano-raw path [path ...]
                        ONT raw reads
  --nano-corr path [path ...]
                        ONT corrected reads
  --subassemblies path [path ...]
                        high-quality contigs input
  -g size, --genome-size size
                        estimated genome size (for example, 5m or 2.6g)
  -o path, --out-dir path
                        Output directory
  -t int, --threads int
                        number of parallel threads [1]
  -i int, --iterations int
                        number of polishing iterations [1]
  -m int, --min-overlap int
                        minimum overlap between reads [auto]
  --asm-coverage int    reduced coverage for initial disjointig assembly [not
                        set]
  --hifi-error float    expected HiFi reads error rate (e.g. 0.01 or 0.001)
                        [0.01]
  --plasmids            rescue short unassembled plasmids
  --meta                metagenome / uneven coverage mode
  --keep-haplotypes     do not collapse alternative haplotypes
  --trestle             enable Trestle [disabled]
  --polish-target path  run polisher on the target sequence
  --resume              resume from the last completed stage
  --resume-from stage_name
                        resume from a custom stage
  --stop-after stage_name
                        stop after the specified stage completed
  --debug               enable debug output
  -v, --version         show program's version number and exit
```

Input reads can be in FASTA or FASTQ format, uncompressed
or compressed with gz. Currently, PacBio (raw, corrected, HiFi)
and ONT reads (raw, corrected) are supported. Expected error rates are
<30% for raw, <3% for corrected, and <1% for HiFi. Note that Flye
was primarily developed to run on raw reads. Additionally, the
--subassemblies option performs a consensus assembly of multiple
sets of high-quality contigs. You may specify multiple
files with reads (separated by spaces). Mixing different read
types is not yet supported. The --meta option enables the mode
for metagenome/uneven coverage assembly.

Genome size estimate is no longer a required option. You
need to provide an estimate if using --asm-coverage option.

To reduce memory consumption for large genome assemblies,
you can use a subset of the longest reads for initial disjointig
assembly by specifying --asm-coverage and --genome-size options. Typically,
40x coverage is enough to produce good disjointigs.

You can run Flye polisher as a standalone tool using
--polish-target option.

## <a name="examples"></a> Examples

You can try Flye assembly on these ready-to-use datasets:

### E. coli P6-C4 PacBio data

The original dataset is available at the 
[PacBio website](https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly).
We coverted the raw `bas.h5` file to the FASTA format for the convenience.

    wget https://zenodo.org/record/1172816/files/E.coli_PacBio_40x.fasta
    flye --pacbio-raw E.coli_PacBio_40x.fasta --out-dir out_pacbio --threads 4

with the `threads` argument being optional (you may adjust it for your environment), 
and `out_pacbio` being the directory where the assembly results will be placed.

### E. coli Oxford Nanopore Technologies data

The dataset was originally released by the 
[Loman lab](http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/).

    wget https://zenodo.org/record/1172816/files/Loman_E.coli_MAP006-1_2D_50x.fasta
    flye --nano-raw Loman_E.coli_MAP006-1_2D_50x.fasta --out-dir out_nano --threads 4


## <a name="inputdata"></a> Supported Input Data

### Oxford Nanopore

We performed our benchmarks with raw ONT reads (R7-R10) with error rate 5-15%.
Due to the biased error pattern, per-nucleotide accuracy is usually lower for 
ONT data than with PacBio data, especially in homopolymer regions. Use `--nano-raw` 
mode for ONT reads. 

### PacBio HiFi

Flye now supports assembly of PacBio HiFi protocol via `--pacbio-hifi` option.
The expected read error is 1% by default. In case the reads are more accurate,
you can adjust `--hifi-error` parameter (for example to 0.001) to potentially
generate more complete assemblies.

### PacBio CLR

Flye was tested on raw PacBio CLR reads (P5C3/P6C4/Sequel) with error rate ~15%.
Use `--pacbio-raw` mode for PacBio CLR.

Note that in PacBio mode, Flye assumes that the input files represent PacBio subreads, 
e.g. adaptors and scraps are removed and multiple passes of the same insertion sequence are separated. 
This is typically handled by PacBio instruments/toolchains, however we saw examples of problemmatic 
raw -> fastq conversions, which resulted into incorrectl subreads. 
In this case, consider using [pbclip](https://github.com/fenderglass/pbclip) to fix your Fasta/q reads.

### Error-corrected reads input

While Flye was designed for assembly of raw reads (and this is the recommended way),
it also supports error-corrected PacBio/ONT reads as input (use the correspondning ```corr``` option).
The parameters are optimized for error rates <3%. If you are getting highly 
fragmented assembly - most likely error rates in your reads are higher. In this case,
consider to assemble using the raw reads instead.

### Consensus of multiple contig sets

```--subassemblies``` input mode generates a consensus of multiple high quality contig assemblies
(such as produced by different short/long read assemblers). The expected error rate
is <1%. You might want to skip the polishing stage with ```--iterations 0``` argument
(however, it might still be helpful to correct small structural errors).


### Input data preparation

Flye works directly with base-called raw reads and does not require any 
prior error correction. Flye automatically detects chimeric reads or reads with low quality ends, 
so you do not need to curate them before the assembly. However, it is always
worth checking for possible contamination in the reads.


## <a name="parameters"></a> Parameter descriptions

### Estimated genome size (optional since 2.8)

No longer required as input. However, it must be used in conjunction with
`--asm-coverage` option.

### Minimum overlap length

This sets a minimum overlap length for two reads to be considered overlapping.
In the latest Flye versions, this parameter is chosen automatically
based on the read length distribution (reads N90) and does not require manual setting.
Typical value is 3k-5k (and down to 1k for datasets with shorter read length).
Intuitively, we want to set this parameter as high as possible, so the
repeat graph is less tangled. However, higher values might lead to assembly gaps.

In some *rare* cases it makes sense to manually increase minimum overlap
for assemblies of big genomes with long reads and high coverage.

### Metagenome mode

Metagenome assembly mode, that is designed for highly non-uniform coverage and
is sensitive to underrepresented sequence at low coverage (as low as 2x).
In some examples of simple metagenomes, we observed that the normal (isolate)
Flye mode assembled more contigious bacterial
consensus sequence, while the metagenome mode was slightly more fragmented, but
revealed strain mixtures. For relatively complex metagenome `--meta` mode
is the recommended way.

### Haplotype mode

By default, Flye (and metaFlye) collapses graph structures caused by
alternative haplotypes (bubbles, superbubbles, roundabouts) to produce
longer consensus contigs. The option `--keep-haplotypes` retains
the alternative paths on the graph, producing less contigouos, but
more detailed assembly.

### Trestle

Trestle is an extra module that resolves simple repeats of
multipicity 2 that were not bridged by reads. Depending on the
datasets, it might resolve a few extra repeats, which is helpful
for small (bacterial genomes). Use `--trestle` option to enable the module.
On large genomes, the contiguity improvements are usually minimal,
but the computation might take a lot of time.

### Reduced contig assembly coverage

Typically, assemblies of large genomes at high coverage require
a hundreds of RAM. For high coverage assemblies, you can reduce memory usage
by using only a subset of longest reads for initial contig extension
stage (usually, the memory bottleneck). The parameter `--asm-coverage`
specifies the target coverage of the longest reads. For a typical assembly, 30x
is enough to produce good initial contigs. Regardless of this parameter,
all reads will be used at the later pipeline stages.

### Number of polishing iterations

Polishing is performed as the final assembly stage. By default, Flye runs one polishing 
iteration. Additional iterations might correct a small number of extra
errors (due to improvements on how reads may align to the corrected assembly). 
If the parameter is set to 0, the polishing is not performed.

### Re-starting from a particular assembly stage

Use `--resume` to resume a previous run of the assembler that may have terminated
prematurely (using the same output directory). 
The assembly will continue from the last previously completed step.

You might also resume from a particular stage with `--resume-from stage_name`,
where `stage_name` is a choice of `assembly, consensus, repeat, trestle, polishing`.
For example, you might supply different sets of reads for different stages.

## <a name="output"></a> Flye output

The main output files are:

* `assembly.fasta` - Final assembly. Contains contigs and possibly scaffolds (see below).
* `assembly_graph.{gfa|gv}` - Final repeat graph. Note that the edge sequences might be
different (shorter) than contig sequences, because contigs might include multiple
graph edges (see below).
* `assembly_info.txt` - Extra information about contigs (such as length or coverage).

Each contig is formed by a single unique graph edge. If possible, unique contigs are
extended with the sequence from flanking unresolved repeats on the graph. Thus,
a contig fully contains the corresponding graph edge (with the same id), but might
be longer then this edge. This is somewhat similar to unitig-contig relation
in OLC assemblers. In a rare case when a repetitive graph edge is not covered by 
the set of "extended" contigs, it will be also output in the assembly file.

Sometimes it is possible to further order contigs into scaffolds based on the 
repeat graph structure. These ordered contigs will be output as a part of scaffold
in the assembly file (with a `scaffold_` prefix).  Since it is hard to give a reliable estimate of the
gap size, those gaps are represented with the default 100 Ns. `assembly_info.txt`
file (below) contains additional information about how scaffolds were formed.

Extra information about contigs/scaffolds is output into the `assembly_info.txt` file.
It is a tab-delimited table with the columns as follows:

* Contig/scaffold id
* Length
* Coverage
* Is circular, (Y)es or (N)o
* Is repetitive, (Y)es or (N)o
* Multiplicity (based on coverage)
* Alternative group
* Graph path (graph path corresponding to this contig/scaffold).

Scaffold gaps are marked with `??` symbols, and `*` symbol denotes a
terminal graph node.

Alternative contigs (representing alternative haplotypes) will have the same
alt. group ID. Primary contigs are marked by `*`

## <a name="graph"></a> Repeat graph

The Flye algorithms are using repeat graph as a core data structure. 
In difference to de Bruijn graphs which require exact k-mer matches,
repeat graphs are built using approximate sequence matches, thus
can tollerate higher noise of SMS reads.

The edges of repeat graph represent genomic sequence, and nodes define
the junctions. All edges are classified into unique and repetitive.
The genome traverses the graph in an unknown way, so as each unique
edge appears exactly once in this traversal. Repeat graphs are useful
for repeat analysis and resolution - which are one of the key 
genome assembly challenges.

<p align="center">
  <img src="graph_example.png" alt="Graph example"/>
</p>

Above is an example of a repeat graph of a bacterial assembly.
Each edge is labeled with its id, length and coverage. Repetitive edges are shown
in color, and unique edges are black. Note that each edge is represented in 
two copies: forward and reverse complement (marked with +/- signs), 
therefore the entire genome is represented in two copies as well. 

In this example, there are two unresolved repeats: (i) a red repeat of 
multiplicity two and length 35k and (ii) a green repeat cluster of multiplicity
three and length 34k - 36k. As the repeats remained unresolved, there are no reads
in the dataset that cover those repeats in full. Five unique edges 
will correspond to five contigs in the final assembly.

Repeat graphs produced by Flye could be visualized using
[AGB](https://github.com/almiheenko/AGB) or [Bandage](https://github.com/rrwick/Bandage).

Repeat graph before repeat resolution could be found in
the `20-repeat/graph_before_rr.gv` file.


## <a name="performance"></a> Flye benchmarks

| Genome                   | Data           | Asm.Size  | NG50     | CPU time  | RAM    |
|--------------------------|----------------|-----------|----------|-----------|--------|
| [E.coli][ecoli]          | PB 50x         | 4.6 Mb    | 4.6 Mb   | 2 h       | 2 Gb   |
| [C.elegans][ce]          | PB 40x         | 106 Mb    | 4.3 Mb   | 100 h     | 31 Gb  |
| [A.thaliana][at]         | PB 75x         | 119 Mb    | 11.9 Mb  | 100 h     | 59 Gb  |
| [D.melanogaster][dm-ont] | ONT 30x        | 136 Mb    | 19.9 Mb  | 130 h     | 33 Gb  |
| [D.melanogaster][dm-pb]  | PB 120x        | 141 Mb    | 18.8 Mb  | 150 h     | 70 Gb  |
| [Human NA12878][na12878] | ONT 35x (rel6) | 2.8 Gb    | 37.9 Mb  | 3100 h    | 394 Gb |
| [Human CHM13 ONT][t2t]   | ONT 120x (rel5)| 2.9 Gb    | 69.4 Mb  | 4000 h    | 450 Gb |
| [Human CHM13 HiFi][t2t]  | PB HiFi 30x    | 3.0 Gb    | 39.8 Mb  | 780 h     | 141 Gb |
| [Human HG002][hg002]     | PB HiFi 30x    | 3.0 Gb    | 33.5 Mb  | 630 h     | 138 Gb |
| [Human CHM1][chm1]       | PB 100x        | 2.8 Gb    | 18.3 Mb  | 2700 h    | 444 Gb |
| [HMP mock][hmp]          | PB meta 7 Gb   | 68 Mb     | 2.6 Mb   | 60 h      | 72 Gb  |
| [Zymo Even][zymo]        | ONT meta 14 Gb | 65 Mb     | 0.7 Mb   | 60 h      | 129 Gb |
| [Zymo Log][zymo]         | ONT meta 16 Gb | 29 Mb     | 0.2 Mb   | 100 h     | 76 Gb  |

[na12878]: https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md
[ce]: https://github.com/PacificBiosciences/DevNet/wiki/C.-elegans-data-set
[at]: https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/
[dm-pb]: https://github.com/PacificBiosciences/DevNet/wiki/Drosophila-sequence-and-assembly
[dm-ont]: https://www.ebi.ac.uk/ena/data/view/SRR6702603
[hg002]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/
[ecoli]: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly
[hmp]: https://github.com/PacificBiosciences/DevNet/wiki/Human_Microbiome_Project_MockB_Shotgun 
[chm1]: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP044331
[t2t]: https://github.com/nanopore-wgs-consortium/CHM13
[zymo]: https://github.com/LomanLab/mockcommunity

The assemblies generated using Flye 2.8 could be downloaded from [Zenodo](https://zenodo.org/record/3965035).
All datasets were run with default parameters for the corresponding read type
with the following exceptions: CHM13 T2T was run with `--min-overlap 10000 --asm-coverage 50`;
CHM1 was run with `--asm-coverage 50`. CHM13 HiFi and HG002 HiFi datasets were run in
`--pacbio-hifi` mode and `--hifi-error 0.003`.

## <a name="algorithm"></a> Algorithm Description

This is a brief description of the Flye algorithm. Please refer to the manuscript
for more detailed information. The draft contig extension is organized as follows:

* K-mer counting / erroneous k-mer pre-filtering
* Solid k-mer selection (k-mers with sufficient frequency, which are unlikely to be erroneous)
* Contig extension. The algorithm starts from a single read and extends it
  with a next overlapping read (overlaps are dynamically detected using the selected
  solid k-mers).

Note that we do not attempt to resolve repeats at this stage, thus
the reconstructed contigs might contain misassemblies. 
Flye then aligns the reads on these draft contigs using minimap2 and
calls a consensus. Afterwards, Flye performs repeat analysis as follows:

* Repeat graph is constructed from the (possibly misassembled) contigs
* In this graph all repeats longer than minimum overlap are collapsed
* The algorithm resolves repeats using the read information and graph structure
* The unbranching paths in the graph are output as contigs

If enabled, after resolving bridged repeats, Trestle module attempts to resolve simple unbridged
repeats (of multiplicity 2) using the heterogeneities between repeat copies.
Finally, Flye performs polishing of the resulting assembly
to correct the remaining errors:

* Alignment of all reads to the current assembly using minimap2
* Partition the alignment into mini-alignments (bubbles)
* Error correction of each bubble using a maximum likelihood approach

The polishing steps could be repeated, which might slightly increase quality for some datasets.
