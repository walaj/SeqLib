[![Build Status](https://travis-ci.org/jwalabroad/SeqLib.svg?branch=master)](https://travis-ci.org/jwalabroad/SeqLib)
[![Coverage Status](https://coveralls.io/repos/github/jwalabroad/SeqLib/badge.svg?branch=master)](https://coveralls.io/github/jwalabroad/SeqLib?branch=master)

C++ interface to HTSlib, BWA-MEM and Fermi

**License:** [Apache2][license]

#NOTE: SeqLib is currently undergoing rapid development changes and some functionality may not be fully implemented.

API Documentation
-----------------
[API Documentation][htmldoc]

Installation
------------

#######
```bash
git clone --recursive https://github.com/jwalabroad/SeqLib.git
cd SeqLib
./configure
make
```
 
I have successfully compiled with GCC-4.8+ on Linux and with Clang on Mac

Description
-----------

SeqLib is a C++ library for querying BAM/SAM/CRAM files, performing 
BWA-MEM operations in memory, and performing sequence assembly. Core operations
in SeqLib are peformed by:
* [HTSlib][htslib]
* [BWA-MEM][BWA] (Apache2 branch)
* [FermiKit][fermi]

SeqLib also has support for storing and manipulating genomic intervals via ``GenomicRegion`` and ``GenomicRegionCollection``. 
It uses an [interval tree][int] (provided by Erik Garrison @ekg) to provide for rapid interval queries.

SeqLib is built to be extendable. See [VariantBam][var] for examples of how to take advantage of C++
class extensions to build off of the SeqLib base functionality. 
 
Memory management
-----------------
SeqLib is built to automatically handle memory management of C code from BWA-MEM and HTSlib by using C++ smart
pointers that handle freeing memory automatically. One of the 
main motivations behind SeqLib is that all access to sequencing reads, BWA, etc should
completely avoid ``malloc`` and ``free``. In SeqLib all the mallocs/frees are handled automatically in the constructors and
destructors.

Note about BamTools, Gamgee and SeqAn
------------------------------
There are overlaps between this project and the [BamTools][BT] project from Derek Barnett, the [Gamgee][gam] 
project from the Broad Institute, and the [SeqAn][seqan] library from Freie Universitat Berlin. These projects 
provide excellent and high quality APIs. SeqLib provides further performance and capabilites for certain classes of 
bioinformatics problems, without attempting to replace these projects.

SeqLib provides some overlapping functionality (eg BAM read/write) but in many cases with improved performance (~2x over BamTools). 
SeqLib further provides in memory access to BWA-MEM, a chromosome aware interval tree and range operations, and to read correction and 
sequence assembly with Fermi. BamTools has more support currently for network access and multi-BAM reading. SeqAn provides 
additional capablities not currently supported in SeqLib, including graph operations and a more expanded suite of multi-sequence alignment
tools (e.g. banded Smith-Waterman). Gamgee provides similar functionality as a C++ interface to HTSlib, but does not incorportate BWA-MEM or Fermi. 
SeqLib is under active development, while Gamgee has been abandoned.

For your particular application, our hope is that SeqLib will provide a comprehensive and powerful envrionment to develop 
bioinformatics tools. Feature requests and comments are welcomed.

Some further SeqLib/BamTools differences
------------------------------
> 1. Sort/index functionality is independently implemented in BamTools. In SeqLib, the Samtools 
 sort and index functions are called directly.
> 2. BamTools stores quality scores and sequences as strings. In SeqLib, the HTSlib ``bam1_t`` format
 is used instead, which uses only 4 bits per base, rather than 8. 
 Conversion to C++ style ``std::string`` is provided with the ``Sequence`` function.
> 3. BamTools provides the ``BamMultiReader`` class for reading multiple BAM files at once, while 
 SeqLib does not currently support this functionality.
> 4. SeqLib contains a built in interface to BWA-MEM for in-memory indexing and querying.
> 5. SeqLib supports reading and writing CRAM files
> 6. BamTools builds with CMake, SeqLib with Autotools.

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
BamRecordVector results;
// hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
bwa.alignSingleSequence("my_seq", querySeq, results, false, 0.9, 10); 

// print results to stdout
for (auto& i : results)
    std::cout << i << std::endl;
## 
```

##### Read a BAM line by line, realign reads with BWA-MEM, write to new BAM
```
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
using SeqLib;

// open the reader BAM/SAM/CRAM
BamReader bw("test.bam");

// open a new interface to BWA-MEM
BWAWrapper bwa;
bwa.retrieveIndex("hg19.fasta");

// open the output BAM
BamWriter writer; // or writer(SeqLib::SAM) or writer(SeqLib::CRAM) 
writer.SetWriteHeader(bwa.HeaderFromIndex());
writer.Open("out.bam");

BamRecord r;
bool rule_passed; // can set rules for what reads to accept. Default is accept all reads. See VariantBam
bool hardclip = false;
float secondary_cutoff = 0.90; // secondary alignments must have score >= 0.9*top_score
int secondary_cap = 10; // max number of secondary alignments to return
while (GetNextRead(r, rule_passed)) {
      BamRecordVector results; // alignment results (can have multiple alignments)
      bwa.alignSingleSequence(r.Sequence(), r.Qname(), results, hardclip, secondary_cutoff, secondary_cap);

      for (auto& i : results)
        writer.writeAlignment(i);
}
```


##### Perform sequence assembly with Fermi directly from a BAM
```

#include "SeqLib/FermiAssembler.h"
using SeqLib;

FermiAssembler f;

// read in data from a BAM
BamReader br("test_data/small.bam");

// retreive sequencing reads (up to 20,000)
BamRecord r;
bool rule;
BamRecordVector brv;
size_t count = 0;
while(br.GetNextRead(r, rule) && count++ < 20000) 
  brv.push_back(r);

// add the reads and error correct them  
f.AddReads(brv);
f.CorrectReads();

// peform the assembly
f.PerformAssembly();

// retrieve the contigs
std::vector<std::string> contigs = f.GetContigs();

// write as a fasta to stdout
for (size_t i = 0; i < contigs.size(); ++i)
    std::cout << ">contig" << i << std::endl << contigs[i] << std::endl;
```

Support
-------
This project is being actively developed and maintained by Jeremiah Wala (jwala@broadinstitute.org)

Attributions
------------
* Jeremiah Wala - Harvard MD-PhD candidate, Bioinformatics and Integrative Genomics, Broad Institute
* Steve Huang - Research Scientist, Broad Institute
* Steve Schumacher - Computational Biologist, Dana Farber Cancer Institute
* Cheng-Zhong Zhang - Research Scientist, Broad Institute
* Marcin Imielinski - Assistant Professor, Cornell University
* Rameen Beroukhim - Assistant Professor, Harvard Medical School

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

[seqan]: https://www.seqan.de

[gam]: https://github.com/broadinstitute/gamgee

[int]: https://github.com/ekg/intervaltree.git

[fermi]: https://github.com/lh3/fermi-lite
