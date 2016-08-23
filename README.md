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

The primary developer for these three projects is Heng Li.

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
bool hardclip = false;
float secondary_cutoff = 0.90; // secondary alignments must have score >= 0.9*top_score
int secondary_cap = 10; // max number of secondary alignments to return
while (GetNextRecord(r)) {
      BamRecordVector results; // alignment results (can have multiple alignments)
      bwa.alignSingleSequence(r.Sequence(), r.Qname(), results, hardclip, secondary_cutoff, secondary_cap);

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
BamReader br("test_data/small.bam");

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

Output:
```
CTATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATCCA                                                                                                    D0EN0ACXX111207:8:2102:10198:18139>>>23:1002968,D0EN0ACXX111207:4:1104:7215:74475>>>23:1003094,
CTATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTAT                    CATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCA                                                D0EN0ACXX111207:4:1202:6670:102657>>>23:1002969,D0UK2ACXX120515:7:2106:5990:22407>>>23:1003056,D0UK2ACXX120515:7:1210:3552:73995>>>23:1003177,
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCA                                                                                            C09DFACXX111207:2:1208:14263:118313>>>23:1002970,C09DFACXX111207:2:2207:16470:21505>>>23:1003098,
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATC                                                                                                      D0EN0ACXX111207:7:2302:16192:143118>>>23:1002971,D0EN0ACXX111207:7:1303:2489:9748>>>23:1003098,
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCA                                                                                            D0EN0ACXX111207:8:1104:7674:93674>>>23:1002972,D0UK2ACXX120515:1:2313:11405:99410>>>23:1003098,
  ATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCAT                                                                                           D0UK2ACXX120515:2:1208:2060:84529>>>23:1002972,C09DFACXX111207:2:1308:4071:198403>>>23:1003099,
   TCTATCTATCTCTTCTTCTGTCCGCTCATGTGTCTGTCCATCTATCTATC                    GTCCATCCATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCCGTCTATCTTATGCATCACAGC                        D0EN0ACXX111207:7:2107:10169:196949>>>23:1002974,D0UK2ACXX120515:2:1213:4675:22715>>>23:1003100,C09DFACXX111207:2:1202:15702:28717>>>23:1003170,
   TCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATT                                                                                          C09DFACXX111207:2:1302:1456:47149>>>23:1002977,D0UK2ACXX120515:1:2210:6216:86040>>>23:1003100,
    CTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTC                                                                                         D0EN0ACXX111207:4:1307:7832:69998>>>23:1002977,D0EN0ACXX111207:8:2201:15817:139174>>>23:1003101,
    CTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTC                                                                                         D0EN0ACXX111207:8:2301:15392:133576>>>23:1002978,C09DFACXX111207:1:1102:1289:63931>>>23:1003101,
      ATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCAT                                                                                       D0UK2ACXX120515:1:1105:14763:99996>>>23:1002980,D0EN0ACXX111207:8:1208:5273:191650>>>23:1003103,
      ATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCAT                                                                                       D0UK2ACXX120515:2:1314:19317:7229>>>23:1002980,D0UK2ACXX120515:6:1312:18321:19662>>>23:1003103,
       TCAATCTCTTCTTCTGTCCGTTAATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTC                                                                                         D0EN0ACXX111207:7:2206:9862:11555>>>23:1002983,D0EN0ACXX111207:4:1208:2280:184516>>>23:1003104,
             TCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATC                                                                                  D0EN0ACXX111207:4:2202:15142:100939>>>23:1002986,D0UK2ACXX120515:6:2103:21074:99867>>>23:1003110,
             TCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCA                                                                                D0EN0ACXX111207:8:1307:1439:161674>>>23:1002986,D0UK2ACXX120515:1:1311:12748:36185>>>23:1003110,
       TCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATC                                                                                      D0UK2ACXX120515:1:2304:14853:22857>>>23:1002986,D0EN0ACXX111207:7:1307:13947:41858>>>23:1003104,
        CTATCTCTTCTTCTGGCCGTTCATGTGTCTGTGCATCTATCTATCCATCTATCTATCATATAACTATCTGTCCATCCATCCATCCATCCA                                                                                                    D0UK2ACXX120515:2:2206:15364:61106>>>23:1002986,D0EN0ACXX111207:8:1102:18033:9414>>>23:1003105,
              CTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCAC                                                                               D0EN0ACXX111207:8:1208:4262:117093>>>23:1002987,D0UK2ACXX120515:6:1113:13735:19005>>>23:1003111,
               TTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACC                                                                              D0UK2ACXX120515:2:2115:2182:97608>>>23:1002987,C09DFACXX111207:2:2206:7329:50468>>>23:1003112,
                 CTTCTGTCCGGTCATGTGTCTGTCCATCTATCTATCCATCTATCTATTATCTAACTA                                                                                                                            D0UK2ACXX120515:2:1313:20357:70197>>>23:1002989,D0EN0ACXX111207:7:1205:5597:147866>>>23:1003114,
                 CTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTA                                                                            D0UK2ACXX120515:6:1208:10024:55268>>>23:1002989,C09DFACXX111207:1:2106:14840:42753>>>23:1003114,
               TTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACC                                                                              D0EN0ACXX111207:4:2103:5351:52080>>>23:1002990,D0UK2ACXX120515:2:2111:15081:99228>>>23:1003112,
               TTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACC                                                                              D0EN0ACXX111207:8:2308:19394:77266>>>23:1002991,D0UK2ACXX120515:1:2311:15980:7405>>>23:1003112,
                    CTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACC                                                                                              D0EN0ACXX111207:7:2303:14877:181503>>>23:1002992,D0UK2ACXX120515:2:2314:11243:85232>>>23:1003117,
                      GTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCAT                                                                       C09DFACXX111207:1:2205:4832:195611>>>23:1002993,D0EN0ACXX111207:7:1201:17225:57893>>>23:1003119,
                      GTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCAT                                                                       D0EN0ACXX111207:8:2203:3295:81886>>>23:1002994,D0UK2ACXX120515:2:2209:9793:61807>>>23:1003119,
                      GTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCAT                                                                       D0UK2ACXX120515:6:2216:15309:35586>>>23:1002994,D0UK2ACXX120515:6:2306:20946:36348>>>23:1003119,
                       TCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATC                                                                      D0UK2ACXX120515:7:1105:7270:23858>>>23:1002994,D0UK2ACXX120515:6:2204:10574:23543>>>23:1003120,
                          GTTCATGTGTCTGTGCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTAT                                                                   D0EN0ACXX111207:7:2106:1668:129298>>>23:1002996,C09DFACXX111207:1:1208:13871:93988>>>23:1003123,
                          GTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTAT                                                                   D0UK2ACXX120515:2:1112:7684:69453>>>23:1002996,D0EN0ACXX111207:8:2106:11477:68328>>>23:1003123,
                           TTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATC                                                                  D0EN0ACXX111207:8:1105:17883:142774>>>23:1002996,C09DFACXX111207:2:1204:6938:118497>>>23:1003124,
C                            CATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAAATATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAA                                                                D0EN0ACXX111207:4:1106:8209:73526>>>23:1002997,D0EN0ACXX111207:4:2202:9848:187980>>>23:1003126,
                             CATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAA                                                                D0UK2ACXX120515:1:1304:11202:38371>>>23:1002997,D0EN0ACXX111207:8:1308:19659:164407>>>23:1003126,
C                             ATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAAT                                                               D0UK2ACXX120515:1:2109:2295:59194>>>23:1002997,D0UK2ACXX120515:6:2207:20734:2444>>>23:1003127,
CTATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTGTC                                                                                                                      D0UK2ACXX120515:7:1201:9029:73248>>>23:1002998,C09DFACXX111207:1:1201:18716:148116>>>23:1003076,
CTA                           ATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAAT                                                               D0UK2ACXX120515:7:1211:9115:89070>>>23:1002999,D0UK2ACXX120515:6:2214:14795:71520>>>23:1003127,
                                GTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCC                                                             D0UK2ACXX120515:7:1102:15726:11108>>>23:1003000,C09DFACXX111207:2:1103:11736:130939>>>23:1003129,
CTC                             GTGTCTGACCATCTATCTATCCATCTTTC                     TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCCG                                            C09DFACXX111207:1:1206:21274:99386>>>23:1003000,C09DFACXX111207:2:2207:19648:153723>>>23:1003129,D0EN0ACXX111207:7:1306:3786:75149>>>23:1003179,
CTCT                              GTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCAT                                                           D0UK2ACXX120515:1:2305:13085:72665>>>23:1003000,D0UK2ACXX120515:6:2202:18237:95154>>>23:1003131,
CTATC                              TCTGTCTATCTATCCATCTATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCAT                                                                       D0EN0ACXX111207:8:1208:8817:16549>>>23:1003001,D0EN0ACXX111207:4:2203:8685:25556>>>23:1003132,
CTATCT                             TCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATC                                                          D0EN0ACXX111207:7:2207:1857:195708>>>23:1003002,D0EN0ACXX111207:8:2103:16744:161185>>>23:1003132,
CTCTTCT                            TCTGTCTATCTATCCATCTATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCAT                                                                       D0UK2ACXX120515:6:1305:17539:97030>>>23:1003003,D0UK2ACXX120515:1:1313:15779:11918>>>23:1003132,
CTATCTA                            TCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTGTCCATCCATC                                                                                                              D0EN0ACXX111207:4:1305:2564:27097>>>23:1003003,D0EN0ACXX111207:7:2307:6657:36792>>>23:1003132,
CTATCTATC                          TCTGTCTATCTATCCATCTATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATACATCCACCCATTCATCCAT                                                                                   D0UK2ACXX120515:6:1214:11566:15357>>>23:1003005,D0EN0ACXX111207:8:2207:1710:15626>>>23:1003132,
CT                                  CTGTCCGTCTATCTATCCATCTATCTATCATCTAA                    CCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCCGTCTATCTTATGCATCACAGCTGGACCATGCAGAAGACA      D0EN0ACXX111207:7:2206:12837:184525>>>23:1003005,C09DFACXX111207:2:1204:15240:111295>>>23:1003133,D0EN0ACXX111207:7:1101:2675:137405>>>23:1003188,
CTATCTATCTA                          TGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATC                                                          D0EN0ACXX111207:8:1301:15362:190820>>>23:1003007,C09DFACXX111207:1:2106:3090:192583>>>23:1003134,
CTATCTATCTAT                         TGTCCATCTATCTATCCATCTATCTATCATCTAAC                    TATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCCGTCTATCTTATGCATCACAGCTGGACCATGCAGAAGACAG     D0UK2ACXX120515:7:1212:11326:2455>>>23:1003008,D0EN0ACXX111207:4:2206:20684:90993>>>23:1003134,D0EN0ACXX111207:7:2206:13050:76135>>>23:1003189,
CT                                    GTCTATCTATCCATCTATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAAT                                                               D0UK2ACXX120515:1:2304:11408:85048>>>23:1003008,D0UK2ACXX120515:6:1101:11721:14638>>>23:1003135,
CTCTTCTTC                             GTCTATCTATCCATCTATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAAT                                                               D0EN0ACXX111207:8:1208:8267:57297>>>23:1003008,D0EN0ACXX111207:7:2106:14449:60792>>>23:1003135,
CTATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATC                    ATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCCGTCTATCTTATGCATCACAGCTGGACCATGCAG            D0EN0ACXX111207:7:2204:17088:122821>>>23:1003009,D0UK2ACXX120515:1:1309:7733:15974>>>23:1003061,D0UK2ACXX120515:1:1307:7793:68870>>>23:1003182,
CTCTTCTTCTGTC                          TCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATC                                                      D0EN0ACXX111207:7:2304:3320:49240>>>23:1003009,D0UK2ACXX120515:2:1201:14912:36409>>>23:1003136,
CTATCTATCTATCT                           CATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCAATCCACCTATCCATCTATCAATCCATC                                                          D0EN0ACXX111207:7:2207:1857:195708>>>23:1003010,C09DFACXX111207:2:2204:18721:94774>>>23:1003138,
CTATCTATCTATCTC                            TCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATC                                                  C09DFACXX111207:2:1102:16071:44068>>>23:1003011,D0UK2ACXX120515:2:2206:1185:3202>>>23:1003140,
CTATCTATCTATCTCTTCTTCTGTCCGT                    CTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCC                                             D0UK2ACXX120515:7:2110:1809:33415>>>23:1003012,C09DFACXX111207:2:2303:20291:10657>>>23:1003145,
CTATCTATCTATCTCTTCTTCTGTCCGTTC                    ATCTATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCAT                                                   D0EN0ACXX111207:8:2203:10017:134410>>>23:1003014,D0UK2ACXX120515:1:2116:7085:14336>>>23:1003147,
CTATCTATCTATCTCTTCTTCTGTCCGTTC                     TCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCCGTC                                          D0UK2ACXX120515:2:1105:17198:38185>>>23:1003014,D0EN0ACXX111207:4:1205:9600:110419>>>23:1003148,
CTATCTATCTATCT                                ATCTATCCATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCAT                                                                                                       C09DFACXX111207:2:1202:10543:32697>>>23:1003014,D0UK2ACXX120515:2:1215:7404:93795>>>23:1003143,
                                              ATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCAT                                               D0EN0ACXX111207:8:2302:20848:61511>>>23:1003015,D0UK2ACXX120515:7:2316:18393:58546>>>23:1003143,
CTATCTATCTATCTCTTCT                            TCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATC                                              D0UK2ACXX120515:6:2315:14105:46897>>>23:1003015,C09DFACXX111207:2:1305:19458:106152>>>23:1003144,
CTATCTATCTATCTCTT                              TCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATC                                              C09DFACXX111207:2:1307:19295:121213>>>23:1003015,D0UK2ACXX120515:7:2310:1519:79721>>>23:1003144,
CTATCTATCTATCTCTTCT                              TATCAATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTAATCCATCCACCTATCCATCTATCA                                                                 D0UK2ACXX120515:7:2204:8522:17890>>>23:1003015,C09DFACXX111207:2:2304:17867:156229>>>23:1003146,

```

Support
-------
This project is being actively developed and maintained by Jeremiah Wala (jwala@broadinstitute.org). 

Attributions
------------
We would like to thank Heng Li (htslib/bwa/fermi), Erik Garrison (interval tree), Christopher Gilbert (aho corasick), 
and Mengyao Zhao (sw alignment), for providing open-source and robust bioinformatics solutions, and Derek Barnett and
the SeqAn team for providing BamTools and SeqAn.

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
