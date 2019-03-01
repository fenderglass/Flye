Flye assembler
==============

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/flye.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/flye)

### Version: 2.4

Flye is a de novo assembler for long and noisy reads, such as
those produced by PacBio and Oxford Nanopore Technologies.
The core read overlapping algorithm uses A-Bruijn graph 
to find the mappings between reads and does not require reads 
to be error-corrected, which resutlts into a faster assembly. 

Flye transforms an initial set of draft contigs into a repeat graph. 
It is somewhat similar to de Bruijn graphs, but can tollerate higher noise level
of long reads. Similarly to short read assemblers, Flye
maps the original reads back on the graph and uses bridging
information to accurately resolve repeats.
The package also includes a polisher module, which produces 
the final assembly of high nucleotide-level quality.


Latest updates
--------------

### Flye 2.4 release (14 Jan 2019)
* Metagenome assembly support fully integrated (`--meta` option)
* New Trestle module for resolving simple unbridged repeats
* New `--plasmids` option that recovers short unassmbled plasmids

### Flye 2.3.7 released (14 Nov 2018)
* Improvements in repeat edges detection
* More precise read mapping - more contigious assemblies for some datasets
* Memory and performance optimizations for high-coverage datasets
* More accurate repeat graphs for complex datasets

### Experimental metagenome version (24 Sep 2018)

Experimental Flye version for metagenome assembly is available 
in the 'flye-meta' branch of this repository. Be aware that
this is still an early 'beta' version, and the code is frequently
updated.


### Flye 2.3.6 released (24 Sep 2018)

* Memory consumption for large genome assemblies reduced by ~30%
* It could be reduced even further by using the new option --asm-coverage,
which specifies a subset of reads for initial contig assembly
* Better repeat graph representation for complex genomes
* Various bugfixes and stability improvements


### Flye 2.3 released (04 Jan 2018)

* ABruijn 2.x branch has been renamed to Flye, highlighting many substantial algorithmic changes
* Stable version of the repeat analysis module
* New command-line syntax (fallback mode with the old syntax is available)
* New --subassemblies mode for generating consensus of multiple assemblies
* Improved preformance and reduced memory footprint (now scales to human genome)
* Corrected reads are now supported
* Extra output with information about the contigs (coverage, multiplicity, graph paths etc.)
* Gzipped Fasta/q support
* Multiple read files support
* Various bugfixes



Manuals
-------

- [Installation instructions](docs/INSTALL.md)
- [Usage](docs/USAGE.md)


Assembly graph
--------------

The Flye algorithms are using repeat graph as a core data structure. 
The edges in this graph represent genomic sequences, and nodes define the
sequence junctions. The genome traverses this graph (in an unknown way) 
so as each unique edge is covered exactly once. The unresolved genomic repeats 
are collapsed into the corresponding edges in the graph
(therefore genome structure remain umbigious).


<p align="center">
  <img src="docs/graph_example.png" alt="Graph example"/>
</p>

Above is an example of a final assembly graph of a bacterial genome.
Each edge is labeled with its id, length and coverage. Repetitive edges are shown
in color, and unique edges are black. The clusters of adjacent repeats are shown with the 
same color. Note that each edge is represented in two copies: forward and
reverse complement (marked with +/- signs), therefore the entire genome is
represented in two copies as well. Sometimes (as in this example), forward and reverse-complement
components are clearly separated, but often they form a single connected component
(in case if the genome contain unresolved inverted repeats).

In this example, there are two unresolved repeats: (i) a red repeat of multiplicity two
and length 35k and (ii) a green repeat cluster of multiplicity three and length 34k - 36k.
As the repeats remained unresolved, there are no reads in the dataset that cover
those repeats in full.


Third-party
-----------

Flye package includes some third-party software:

* [libcuckoo](http://github.com/efficient/libcuckoo)
* [intervaltree](https://github.com/ekg/intervaltree)
* [lemon](http://lemon.cs.elte.hu/trac/lemon)
* [minimap2](https://github.com/lh3/minimap2)


License
-------

Flye is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Flye is developed in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)

Code contributions:

* Repeat graph and current package maintaining: Mikhail Kolmogorov
* Trestle module and original polisher code: Jeffrey Yuan
* Original contig extention code: Yu Lin
* Short plasmids recovrey module: Evgeny Polevikov


Publications
------------
Mikhail Kolmogorov, Jeffrey Yuan, Yu Lin and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using Repeat Graphs", bioRxiv, 2018

Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Mark Chaisson and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs", PNAS, 2016


Contacts
--------
Please report any problems to the github issue tracker (at http://github.com/fenderglass/Flye).
If possible, please include "flye.log" file from the output directory
for faster feedback. Alternatively, you can write directly to fenderglass@gmail.com.
