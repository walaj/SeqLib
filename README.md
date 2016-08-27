[![Build Status](https://travis-ci.org/jwalabroad/SeqLib.svg?branch=master)](https://travis-ci.org/jwalabroad/SeqLib)
[![Coverage Status](https://coveralls.io/repos/github/jwalabroad/SeqLib/badge.svg?branch=master)](https://coveralls.io/github/jwalabroad/SeqLib?branch=master)

C++ interface to HTSlib, BWA-MEM and Fermi

**License:** 
-- This branch incorporates source code from Blat, which is copyrighted by Jim Kent and with a [license][blatlicense] for academic, personal and non-profit use. To
use this branch for commercial use, see the Blat license in SeqLib/blat/README.
-- For an Apache2 version without Blat integration, see the master branch.

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

Integrating into build system
-----------------------------

After building, you will need to add the relevant header directories:
```bash
SEQ=<path_to_seqlib>
C_INCLUDE_PATH=$C_INCLUDE_PATH:$SEQ:$SEQ/htslib:$SEQ/src:$SEQ/fermi-lite
```

And need to link the following libraries
```bash
SEQ=<path_to_seqlib>
LDADD="$LDADD -L$SEQ/src/libseqlib.a -L$SEQ/fermi-lite/libfml.a -L$SEQ/libbwa.a -L$SEQ/htslib/libhts.a"
```

Description
-----------

SeqLib is a C++ library for querying BAM/SAM/CRAM files, performing 
BWA-MEM operations in memory, and performing sequence assembly. Core operations
in SeqLib are peformed by:
* [HTSlib][htslib]
* [BWA-MEM][BWA] (Apache2 branch)
* [Blat][blatlink] 
* [FermiKit][fermi]

The primary developer for the HTSlib, BWA-MEM and FermiKit projects is Heng Li. Blat is developed by Jim Kent.

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

Example usages
--------------
##### Targeted re-alignment of reads to a given region with BWA-MEM
```
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
using SeqLib;
RefGenome ref;
ref.LoadIndex("hg19.fasta");

// get sequence at given locus
std::string seq = ref.queryRegion("1", 1000000,1001000);

// Make an in-memory BWA-MEM index of region
BWAWrapper bwa;
UnalignedSequenceVector usv = {{"chr_reg1", seq}};
bwa.ConstructIndex(usv);

// align an example string with BWA-MEM
std::string querySeq = "CAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAG";
BamRecordVector results;
// hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
bwa.AlignSequence("my_seq", querySeq, results, false, 0.9, 10); 

// print results to stdout
for (auto& i : results)
    std::cout << i << std::endl;
```

##### Read a BAM line by line, realign reads with BWA-MEM, write to new BAM
```
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
using SeqLib;

// open the reader BAM/SAM/CRAM
BamReader bw;
bw.Open("test.bam");

// open a new interface to BWA-MEM
BWAWrapper bwa;
bwa.LoadIndex("hg19.fasta");

// open the output BAM
BamWriter writer; // or writer(SeqLib::SAM) or writer(SeqLib::CRAM) 
writer.SetWriteHeader(bwa.HeaderFromIndex());
writer.Open("out.bam");
writer.WriteHeader();

BamRecord r;
bool hardclip = false;
float secondary_cutoff = 0.90; // secondary alignments must have score >= 0.9*top_score
int secondary_cap = 10; // max number of secondary alignments to return
while (GetNextRecord(r)) {
      BamRecordVector results; // alignment results (can have multiple alignments)
      bwa.AlignSequence(r.Sequence(), r.Qname(), results, hardclip, secondary_cutoff, secondary_cap);

      for (auto& i : results)
        writer.WriteRecord(i);
}
```


##### Perform sequence assembly with Fermi directly from a BAM
```

#include "SeqLib/FermiAssembler.h"
using SeqLib;

FermiAssembler f;

// read in data from a BAM
BamReader br;
br.Open("test_data/small.bam");

// retreive sequencing reads (up to 20,000)
BamRecord r;
BamRecordVector brv;
size_t count = 0;
while(br.GetNextRead(r) && count++ < 20000) 
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

##### Plot a collection of gapped alignments
```
using SeqLib;
BamReader r;
r.Open("test_data/small.bam");

GenomicRegion gr("X:1,002,942-1,003,294", r.Header());
r.SetRegion(gr);

SeqPlot s;
s.SetView(gr);

BamRecord rec;
BamRecordVector brv;
while(r.GetNextRecord(rec))
  if (!rec.CountNBases() && rec.MappedFlag())
    brv.push_back(rec);
s.SetPadding(20);

std::cout << s.PlotAlignmentRecords(brv);
```

Trimmed output from above (NA12878):
```
CTATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATCCA
CTATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTAT                    CATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCA
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCA              
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATC                        
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCA              
  ATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCAT             
   TCTATCTATCTCTTCTTCTGTCCGCTCATGTGTCTGTCCATCTATCTATC                    GTCCATCCATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCCGTCTATCTTATGCATCACAGC
   TCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATT                                                                  
    CTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTC                                                                 
    CTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTC                                                                 
      ATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCAT                                                               
      ATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCAT                                                               
```

##### Read simultaneously from a BAM, CRAM and SAM file and send to stdout
```
using SeqLib;
#include "SeqLib/BamReader.h"
BamReader r;
  
// read from multiple streams coordinate-sorted order
r.Open("test_data/small.bam");
r.Open("test_data/small.cram");
r.Open("test_data/small.sam");

BamWriter w(SeqLib::SAM); // set uncompressed output
w.Open("-");              // write to stdout
w.SetHeader(r.Header());  // specify the header
w.WriteHeader();          // write out the header

BamRecord rec;
while(r.GetNextRecord(rec)) 
  w.WriteRecord(rec);
w.Close();               // Optional. Will close on destruction
```

Support
-------
This project is being actively developed and maintained by Jeremiah Wala (jwala@broadinstitute.org). 

Attributions
------------
We would like to thank Heng Li (htslib/bwa/fermi), Erik Garrison (interval tree), Christopher Gilbert (aho corasick), 
and Mengyao Zhao (sw alignment), for providing open-source and robust bioinformatics solutions.

Development, support, guidance, testing:
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

[blatlicense]: https://github.com/jwalabroad/SeqLib/blob/old_license/blat/README

[blatlink]: http://genome.cshlp.org/content/12/4/656.full
