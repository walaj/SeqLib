[![Build Status](https://travis-ci.org/jwalabroad/SeqLib.svg?branch=master)](https://travis-ci.org/jwalabroad/SeqLib)
[![Coverage Status](https://coveralls.io/repos/github/jwalabroad/SeqLib/badge.svg?branch=master)](https://coveralls.io/github/jwalabroad/SeqLib?branch=master)

C++ interface to HTSlib, BWA-MEM and Fermi-kit

**License:** [Apache2][license]

API Documentation
-----------------
[API Documentation][htmldoc]

Installation
------------

#######
```bash
git clone https://github.com/jwalabroad/SeqLib.git
cd SeqLib
./configure
make
```
 
I have successfully compiled with GCC-4.8+ on Linux and with Clang on Mac

Description
-----------

SeqLib is a C++ package for querying BAM and SAM files, performing 
BWA-MEM and BLAT operations in memory, and performing advanced filtering of 
reads using a hierarchy of rules. Currently, SeqLib wraps the following projects:
* [HTSlib][htslib]
* [BWA-MEM][BWA] (Apache2 branch)
* [FermiKit][fermi] (Apache2 branch)

SeqLib also has support for storing and manipulating genomic intervals via ``GenomicRegion`` and ``GenomicRegionCollection``. 
It uses an [interval tree][int] (provided by Erik Garrison @ekg) to provide for rapid interval queries.

SeqLib is built to be extendable. See [Variant Bam][var] for examples of how to take advantage of C++
class extensions to build off of the SeqLib base functionality. 
 
Memory management
-----------------
SeqLib is built to automatically handle memory management of C code from BWA-MEM and htslib by using C++ smart
pointers that handle freeing memory automatically. One of the 
main motivations behind SeqLib is that all access to sequencing reads, BWA, etc should
completely avoid ``malloc`` and ``free``. In SeqLib, the speed and compression of HTSlib
is available, but all the mallocs/frees are handled automatically in the constructors and
destructors.

Note about BamTools and Gamgee
------------------------------
There are overlaps between this project and the [BamTools][BT] project from Derek Barnett, and the [Gamgee][gam] 
project from the Broad Institute. In short, BamTools has been more widely used and tested, but is relatively slow compared with SeqLib (~2x).
Gamgee provides similar functionality as a C++ interface to HTSlib, but does not incorportate BWA-MEM or BLAT. SeqLib is under active development, while Gamgee
has been abandoned.

SeqLib/BamTools differences
------------------------------
> 1. Sort/index functionality is independently implemented in BamTools. In SeqLib, the Samtools 
 sort and index functions are called directly.
> 2. BamTools stores quality scores and sequences as strings. In SeqLib, the HTSlib ``bam1_t`` format
 is used instead, which uses only 4 bits per base, rather than 8. 
 Conversion to C++ style ``std::string`` is provided with the ``Sequence`` function.
> 3. BamTools provides the ``BamMultiReader`` class for reading multiple BAM files at once, while 
 SeqLib does not currently support this functionality.
> 4. SeqLib contains a built in interface to BWA-MEM for in-memory indexing and querying.
> 5. SeqLib contains a beta wrapper around BLAT.
> 6. SeqLib supports reading and writing CRAM files
> 7. BamTools has been widely used in a number of applications, and is thus substantially more tested.
> 8. SeqLib is faster at reading/writing BAM files by about 2x.
> 9. BamTools builds with CMake, SeqLib with Autotools.

Example usages
--------------
##### Targeted re-alignment of reads to a given region with BWA-MEM
```
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
using SeqLib;
RefGgenome ref("hg19.fasta");

## get sequence at given locus
std::string seq = ref.queryRegion("1", 1000000,1001000);

## Make an in-memory BWA-MEM index of region
BWAWrapper bwa;
USeqVector usv = {{"chr_reg1", seq}};
bwa.constructIndex(usv);

## align an example string with BWA-MEM
std::string querySeq = "CAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAG";
BamReadVector results;
// hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
bwa.alignSingleSequence("my_seq", querySeq, results, false, 0.9, 10); 

// print results to stdout
for (auto& i : results)
    std::cout << i << std::endl;
## 
```

##### Read a BAM line by line, realign reads with BWA-MEM, write to new BAM
```
#include "SeqLib/BamWalker.h"
#include "SeqLib/BWAWrapper.h"
using SeqLib;

// open the reader BAM
BamWalker bw("test.bam");

// open a new interface to BWA-MEM
BWAWrapper bwa;
bwa.retrieveIndex("hg19.fasta");

// open the output bam
BamWalker writer;
write.SetWriteHeader(bwa.HeaderFromIndex());
write.OpenWriteBam("out.bam");

BamRead r;
bool rule_passed; // can set rules for what reads to accept. Default is accept all reads. See VariantBam
bool hardclip = false;
float secondary_cutoff = 0.90; // secondary alignments must have score >= 0.9*top_score
int secondary_cap = 10; // max number of secondary alignments to return
while (GetNextRead(r, rule_passed)) {
      BamReadVector results; // alignment results (can have multiple alignments)
      bwa.alignSingleSequence(r.Sequence(), r.Qname(), results, hardclip, secondary_cutoff, secondary_cap);

      for (auto& i : results)
        writeAlignment(i);
}
```


Support
-------
This project is being actively developed and maintained by Jeremiah Wala (jwala@broadinstitute.org)

Attributions
------------
* Jeremiah Wala - Harvard MD-PhD candidate, Bioinformatics and Integrative Genomics, Broad Institute
* Steve Huang - Scientist, Broad Institute
* Steve Schumacher - Research Scientist, Dana Farber Cancer Institute
* Cheng-Zhong Zhang - Research Scientist, Broad Institute
* Marcin Imielinski - Assistant Professor, Cornell University

[htslib]: https://github.com/samtools/htslib.git

[SGA]: https://github.com/jts/sga

[BLAT]: https://genome.ucsc.edu/cgi-bin/hgBlat?command=start

[BWA]: https://github.com/lh3/bwa

[license]: https://github.com/jwalabroad/SeqLib/blob/new_license/LICENSE

[BamTools]: https://raw.githubusercontent.com/wiki/pezmaster31/bamtools/Tutorial_Toolkit_BamTools-1.0.pdf

[API]: http://pezmaster31.github.io/bamtools/annotated.html

[htmldoc]: http://jwalabroad.github.io/SeqLib/doxygen

[var]: https://github.com/jwalabroad/VariantBam

[BT]: https://github.com/pezmaster31/bamtools

[gam]: https://github.com/broadinstitute/gamgee

[int]: https://github.com/ekg/intervaltree.git

[fermi]: https://github.com/lh3/fermi-lite
