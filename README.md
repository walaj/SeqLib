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
git clone --recursive https://github.com/boostorg/boost.git
cd boost
./bootstrap.sh --with-libraries=regex
./b2


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
* [SGA][SGA]
* [BLAT][BLAT]

SnowTools is built to be extendable. See [Variant Bam][var] for examples of how to take advantage of C++
class extensions to build off of the SnowTools base functionality. 

As noted above, there are many overlaps between this project at the [BamTools][BT] project from Derek Barnett.
SnowTools is provided here in the case that it may be more suited to your individual needs than BamTools, or vice-versa. To
aid developers in deciding which package is right for them, I have put together a small list of 
similarities and differences between these two packages
 
Example usages
--------------
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

SnowTools/BamTools similarity
-----------------------------

> 1. Provide read/write access to BAM files
> 2. Classes containing individual reads and ways to interact with them (e.g. edit tags)
> 3. Primarily documented through Doxygen docs

SnowTools/BamTools differences
------------------------------
SnowTools/BamTools differences
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

Support
-------
This code is being developed and maintained by Jeremiah Wala (jwala@broadinstitute.org)

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
