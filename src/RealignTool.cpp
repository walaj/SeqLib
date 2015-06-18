#include "RealignTool.h"

#include <string>
#include <iostream>

#include "SnowTools/AlignedContig.h"
#include "SnowTools/BamWalker.h"

#define _DIVBWT 1
#define FILE_DUMP 1

namespace opt {

  std::string refgenome = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  std::string bamin = "/seq/picard_aggregation/G14856/TCGA-05-4432-01A-01D-1931-08/v2/TCGA-05-4432-01A-01D-1931-08.bam";
}

void runRealign(int argc, char** argv) 
{

  // make a BAMWalker
  SnowTools::BamWalker bw(opt::bamin);
  bw.setStdout();
  std::cout << bw << std::endl;
  
  SnowTools::USeqVector v = {
    {"TP53", "GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCTTCTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGGTAAGCTCCTGACTGAACTTGATGAGTCCTCTCTGA"},
    {"MYC", "GACCCCCGAGCTGTGCTGCTCGCGGCCGCCACCGCCGGGCCCCGGCCGTCCCTGGCTCCCCTCCTGCCTCGAGAAGGGCAGGGCTTCTCAGAGGCTTGGCGGGAAAAAGAACGGAGGGAGGGATCGCGCTGAGTATAAAAGCCGGTTTTCGGGGCTTTATCTAACTCGCTGTAGTAATTCCAGCGAGAGGCAGAGGGAGCGAGCGGGCGG"}
  };

  SnowTools::BWAWrapper w;
  //w.constructIndex(v);
  w.retrieveIndex(opt::refgenome);

  //std::string test_seq = v[0].seq + v[1].seq;
  std::string test_seq = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTA";
  
  std::string name = "TP53_MYC";

  SnowTools::BamReadVector brv;
  w.alignSingleSequence(test_seq, name, brv, false);
  return;

  SnowTools::AlignedContig ac(brv);
  
  for (auto& i : brv)
    bw.WriteAlignment(i);
  
  // get 
  SnowTools::GenomicRegionVector regs;
  regs.push_back(SnowTools::GenomicRegion("17:7489708-7591708", bw.header()));
  regs.push_back(SnowTools::GenomicRegion("17:7489708-7591708", bw.header()));
  bw.setBamWalkerRegions(regs);
  SnowTools::BamRead b;
  bool rule;
  SnowTools::BamReadVector read_store;
  while (bw.GetNextRead(b, rule)) 
    {
      b.AddZTag("SR", "t" + std::to_string(b.AlignmentFlag()) + "_" + b.Qname());
      read_store.push_back(b);
    }

  ac.alignReads(read_store);

  std::cerr << ac << std::endl;

  
}
