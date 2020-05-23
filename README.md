Flye assembler
==============

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/flye.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/flye)

### Version: 2.7.1

Flye is a de novo assembler for single molecule sequencing reads,
such as those produced by PacBio and Oxford Nanopore Technologies.
It is designed for a wide range of datasets, from small bacterial projects
to large mammalian-scale assemblies. The package represents a complete
pipeline: it takes raw PacBio / ONT reads as input and outputs polished contigs.
Flye also has a special mode for metagenome assembly.

Manuals
-------

- [Installation instructions](docs/INSTALL.md)
- [Usage](docs/USAGE.md)
- [FAQ](docs/FAQ.md)

Latest updates
--------------

### Flye 2.7.1. release (24 Apr 2020)
* Fixes very long GFA generation time for some large assemblies (no other changes)

### Flye 2.7 release (03 Mar 2020)
* Better assemblies of real (and comlpex) metagenomes
* New option to retain alternative haplotypes, rather than collapsing them (`--keep-haplotypes`)
* PacBio HiFi mode
* Using Bam instead of Sam to reduce storage requirements and IO load
* Improved human assemblies
* Annotation of alternative contigs
* Better polishing quality for the newest ONT datasets
* Trestle module is disabled by default (use `--trestle` to enable)
* Many big fixes and improvements

### Flye 2.6 release (19 Sep 2019)
* This release introduces Python 3 support (no other changes)

### Flye 2.5 release (25 Jul 2019)
* Better ONT polishing for the latest basecallers (Guppy/flipflop)
* Improved consensus quality of repetitive regions
* More contiguous assemblies of real metagenomes
* Improvements for human genome assemblies
* Various bugfixes and performance optimizations
* Also check the new [FAQ section](docs/FAQ.md)


Repeat graph
------------

Flye is using repeat graph as a core data structure. 
In difference to de Bruijn graphs (which require exact k-mer matches),
repeat graphs are built using approximate sequence matches, and
can tolerate higher noise of SMS reads.

The edges of repeat graph represent genomic sequence, and nodes define
the junctions. Each edges is classified into unique or repetitive.
The genome traverses the graph (in an unknown way), so as each unique
edge appears exactly once in this traversal. Repeat graphs reveal the
repeat structure of the genome, which helps to reconstruct an optimal assembly.


<p align="center">
  <img src="docs/graph_example.png" alt="Graph example"/>
</p>

Above is an example of the repeat graph of a bacterial assembly.
Each edge is labeled with its id, length and coverage. Repetitive edges are shown
in color, and unique edges are black. Note that each edge is represented in 
two copies: forward and reverse complement (marked with +/- signs), 
therefore the entire genome is represented in two copies. This is necessary
because the orientation of input reads is unknown.

In this example, there are two unresolved repeats: (i) a red repeat of 
multiplicity two and length 35k and (ii) a green repeat cluster of multiplicity
three and length 34k - 36k. As the repeats remained unresolved, there are no reads
in the dataset that cover those repeats in full. Five unique edges 
will correspond to five contigs in the final assembly.

Repeat graphs produced by Flye could be visualized using
[AGB](https://github.com/almiheenko/AGB) or [Bandage](https://github.com/rrwick/Bandage).


Flye benchmarks
---------------

| Genome                   | Data           | Asm.Size  | NG50     | CPU time  | RAM    |
|--------------------------|----------------|-----------|----------|-----------|--------|
| [E.coli][ecoli]          | PB 50x         | 4.6 Mb    | 4.6 Mb   | 2 h       | 2 Gb   |
| [C.elegans][ce]          | PB 40x         | 102 Mb    | 3.6 Mb   | 100 h     | 31 Gb  |
| [A.thaliana][at]         | PB 75x         | 120 Mb    | 9.5 Mb   | 100 h     | 46 Gb  |
| [D.melanogaster][dm-ont] | ONT 30x        | 139 Mb    | 10.6 Mb  | 130 h     | 31 Gb  |
| [D.melanogaster][dm-pb]  | PB 120x        | 142 Mb    | 18.8 Mb  | 150 h     | 75 Gb  |
| [Human NA12878][na12878] | ONT 35x (rel6) | 2.9 Gb    | 33.2 Mb  | 2500 h    | 714 Gb |
| [Human CHM13 T2T][t2t]   | ONT 120x (rel3)| 2.9 Gb    | 75.1 Mb  | 5000 h    | 871 Gb |
| [Human HG002][hg002]     | PB CCS 30x     | 2.9 Gb    | 27.5 Mb  | 1400 h    | 272 Gb |
| [Human CHM1][chm1]       | PB 100x        | 2.8 Gb    | 21.5 Mb  | 2700 h    | 676 Gb |
| [HMP mock][hmp]          | PB meta 7 Gb   | 66 Mb     | 2.6 Mb   | 60 h      | 72 Gb  |
| [Zymo Even][zymo]        | ONT meta 14 Gb | 64 Mb     | 0.6 Mb   | 60 h      | 129 Gb |
| [Zymo Log][zymo]         | ONT meta 16 Gb | 23 Mb     | 1.3 Mb   | 100 h     | 76 Gb  |

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

The assemblies generated using Flye 2.7 could be downloaded from [Zenodo](https://zenodo.org/record/3694400).
All datasets were run with default parameters for the corresponding read type
with the following exceptions: CHM13 T2T was run with `--min-overlap 10000 --asm-coverage 50`;
CHM1 was run with `--asm-coverage 40`.

Third-party
-----------

Flye package includes some third-party software:

* [libcuckoo](http://github.com/efficient/libcuckoo)
* [intervaltree](https://github.com/ekg/intervaltree)
* [lemon](http://lemon.cs.elte.hu/trac/lemon)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](https://https://github.com/samtools/samtools)


License
-------

Flye is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Flye is developed in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)

Code contributions:

* Repeat graph and current package maintaining: Mikhail Kolmogorov
* Trestle module and original polisher code: Jeffrey Yuan
* Original contig extension code: Yu Lin
* Short plasmids recovery module: Evgeny Polevikov


Publications
------------
Mikhail Kolmogorov, Jeffrey Yuan, Yu Lin and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using Repeat Graphs", Nature Biotechnology, 2019
[doi:10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8)

Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Mark Chaisson and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs", PNAS, 2016
[doi:10.1073/pnas.1604560113](https://www.doi.org/10.1073/pnas.1604560113)


How to get help
---------------
A preferred way report any problems or ask questions about Flye is the 
[issue tracker](https://github.com/fenderglass/Flye/issues). 
Before posting an issue/question, consider to look through the FAQ
and existing issues (opened and closed) - it is possble that your question
has already been answered.

If you reporting a problem, please include the `flye.log` file and provide
details about your dataset.

In case you prefer personal communication, please contact Mikhail at fenderglass@gmail.com.
