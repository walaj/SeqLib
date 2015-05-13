#include "RealignTool.h"

#include <string>
#include <iostream>

#include "SnowTools/BWAWrapper.h"

#define _DIVBWT 1
#define FILE_DUMP 1

namespace opt {

  std::string refgenome = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
}

void runRealign(int argc, char** argv) 
{

  SnowTools::USeqVector v = {
    {"test1", "TGCCACTGTACTCCAGCTTGGACAACAGACTGAGACCCTGTCTCGAAAGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAGAAAGAAAGAAAGAAAGAAAAGG"},
    {"test2", "TGCCACTGTACTCCAGCTTGGACAACAGACTGAGACCCTGTCTCGAAAGAAGGAAGGAAGGAAGGAAGGAA"},
    {"test3", "ACTCCAGCTTGGACAACAACTGAGACCCTGTCTCGAAAGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAGAAAGAAAGAAAGAA"}
  };

  SnowTools::BWAWrapper w;
  w.constructIndex(v);
  w.writeIndexToFiles("myidx.fa");
  //w.retrieveIndex("ref.fa");//opt::refgenome);

  std::string test_seq = "ACTCCAGCTTGGACAACAACTGAGACCCTGTCTCGAAAGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAGAAAGAAAGAAAGAA";
  w.alignSingleSequence(test_seq);
  
  std::cout << "DONE" << std::endl;
}
