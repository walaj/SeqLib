[![Build Status](https://travis-ci.org/jwalabroad/SnowTools.svg?branch=master)](https://travis-ci.org/jwalabroad/SnowTools)
[![Coverage Status](https://coveralls.io/repos/jwalabroad/SnowTools/badge.svg?branch=master&service=github)](https://coveralls.io/github/jwalabroad/SnowTools?branch=master)

<div style="text-align:center"><img src="https://raw.githubusercontent.com/jwalabroad/SnowTools/master/figs/labels_st.png" width="250"></div>

C++ htslib/bwa-mem interface and command line tools for interrogating BAM and SAM files.

**License:** [GNU GPLv3][license]

API Documentation
-----------------
[API Documentation][htmldoc]

Installation
------------

#######
```bash
### if on Broad Institute servers, add GCC-4.9
reuse -q GCC-4.9

############## DOWNLOAD AND INSTALL BOOST ###############
############## (only if not already installed) ##########
wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz
tar -xvzf boost_1_61_0.tar.gz
## we only user header-only libraries, so no compiling of Boost is needed

############### DOWNLOAD SNOWTOOLS ############### 
git clone https://github.com/jwalabroad/SnowTools.git
cd SnowTools

############### COMPILE AND INSTALL ###############
./configure --with-boost=<path_to_boost>
make
```
 
I have successfully compiled with GCC-4.8+ on Linux.

Description
-----------

SnowTools is a C++ package for querying BAM and SAM files, performing 
BWA-MEM and BLAT operations in memory, and performing advanced filtering of 
reads using a hierarchy of rules. Currently, SnowTools wraps the following projects:
* [HTSlib][htslib]
* [BWA-MEM][BWA]
* [BLAT][BLAT]

SnowTools also has support for storing and manipulating genomic intervals via ``GenomicRegion`` and ``GenomicRegionCollection``. 
It uses an [interval tree][int] (provided by Erik Garrison @ekg) to provide for rapid interval queries.

SnowTools is built to be extendable. See [Variant Bam][var] for examples of how to take advantage of C++
class extensions to build off of the SnowTools base functionality. 
 
Memory management
-----------------
One of the greater challenges in using C code like BWA-MEM and htslib 
as an API is in handling memory management. C++ makes this more palatable with smart
pointers that handle freeing memory automatically. One of the 
main motivations behind SnowTools is that all access to sequencing reads, BWA, etc should
completely avoid ``malloc`` and ``free``. In SnowTools, the speed and compression of HTSlib
is available, but all the mallocs/frees are handled for you automatically in the constructors and
destructors.

Note about BamTools and Gamgee
------------------------------
There are many overlaps between this project and the [BamTools][BT] project from Derek Barnett, and the [Gamgee][gam] 
project also from the Broad Institute. These are excellent projects, and I am providing SnowTools here 
in the case that it may be more suited to your individual needs than BamTools or Gamgee. 

In short, BamTools is the mostly widely used and tested program from BAM/SAM reading and writing in C++, but 
is ~2x slower than htslib and has a larger memory footprint. Gamgee provides BAM/SAM/CRAM reading/writing, 
and uses HTSlib to maximize efficiency, plus has smart pointer memory management. It is a more mature project, but does not have an interface to BWA-MEM or BLAT.
SnowTools is the least mature of these projects, but unites in one eco system HTSlib, BWA-MEM, BLAT, and soon SGA (String Graph Assembler).

SnowTools/BamTools differences
------------------------------
> 1. Sort/index functionality is independently implemented in BamTools. In SnowTools, the Samtools 
 sort and index functions are called directly.
> 2. BamTools stores quality scores and sequences as strings. In SnowTools, the HTSlib native bam1_t format
 is used instead. This format has a lower memory footprint by using only 4 bits per base, rather than 8. 
 Conversion to C++ style std::string is provided as a function and can be done on the fly.
> 3. BamTools provides the BamMultiReader class for reading multiple BAM files at once, while 
 SnowTools does not currently support this functionality.
> 4. SnowTools contains a built in interface to BWA-MEM for in-memory indexing and querying.
> 5. SnowTools contains a beta wrapper around BLAT.
> 6. SnowTools supports reading and writing CRAM files
> 7. BamTools has been widely used in a number of applications, and is thus substantially more tested.
> 8. SnowTools is faster at reading/writing BAM files by about 2x.
> 9. BamTools builds with CMake, SnowTools with Autotools.

Example usages
--------------
##### Targeted re-alignment of reads to a given region with BWA-MEM
```
#include "SnowTools/RefGenome.h"
#include "SnowTools/BWAWrapper.h"
using SnowTools;
RefGgenome ref("hg19.fasta");

## get sequence at given locus
std::string seq = ref.queryRegion("1", 1000000,1001000);

## Make an in-memory BWA-MEM index of region
BWAWrapper bwa;
USeqVector usv = {{"chr_reg1", seq}};
bwa.constructIndex(usv);

## query the string
std::string querySeq = "CAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAG";
BamReadVector results;
// hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
bwa.alignSingleSequence("my_seq", querySeq, results, false, 0.9, 10); 

// print results to stdout
for (auto& i : results)
    std::cout << i << std::endl;
## 
```

##### Read a BAM line by line, realign reads with BWA-MEM
```
#include "SnowTools/BamWalker.h"
#include "SnowTools/BWAWrapper.h"
using SnowTools;

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
bool rule; // can set rules for what reads to accept. Default is all, so rule always true.
bool hardclip = false;
float secondary_cutoff = 0.90; // secondary alignments must have score >= 0.9*top_score
int secondary_cap = 10; // max number of secondary alignments to return
while (GetNextRead(r, rule)) {
      BamReadVector results; // alignment results (can have multiple alignments)
      bwa.alignSingleSequence(r.Sequence(), r.Qname(), results, hardclip, secondary_cutoff, secondary_cap);

      for (auto& i : results)
        write.Alignment(i);
}
```


Support
-------
This code is being actively developed and maintained by Jeremiah Wala (jwala@broadinstitute.org)

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

[license]: https://github.com/broadinstitute/variant-bam/blob/master/LICENSE

[BamTools]: https://raw.githubusercontent.com/wiki/pezmaster31/bamtools/Tutorial_Toolkit_BamTools-1.0.pdf

[API]: http://pezmaster31.github.io/bamtools/annotated.html

[htmldoc]: http://jwalabroad.github.io/SnowTools/doxygen

[var]: https://github.com/jwalabroad/VariantBam

[BT]: https://github.com/pezmaster31/bamtools

[gam]: https://github.com/broadinstitute/gamgee

[int]: https://github.com/ekg/intervaltree.git
