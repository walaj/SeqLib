#include <string>
#include <getopt.h>
#include <iostream>

#include "SnowTools/BLATWrapper.h"

int main(int argc, char** argv) {

  SnowTools::BLATWrapper bw;

  std::string ooc = "/xchip/gistic/Jeremiah/blat/11.ooc";

  std::cerr << "...loading BLAT index" << std::endl;
  bw.loadIndex("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta", ooc);
  
  std::cerr << "...querying file" << std::endl;
  SnowTools::BamReadVector brv, brv2;
  std::string seq = "ATGGTGGTAGTTAGGGTCATGGTGGTAGTTAGGGTCACGGTGGTGGTTAGGGTCACGGTGGTGGTTAGAGTCACAGGGTAGAACCCTTGTGGTGGGATTTGTGCCCTTTATAGGATGAGAGGATGAGACACAAGAGAGGTTGTGCTGCGCCTGTGCTCTCTGCTCCACATGAGAACATGGTGAGCATGAGGCCGCCAGCAAGCAAGGAGATACCCCGCCCTGCAGGTTCCGTCATCCTGACTCCAGCCTCGGAAACATGAGAAAGTCAATGCCTGTCACTTAAGCCGCCCAGTCTGTGGTATTTTGCTGTGGTGGCTGAGCCGACGGGGACAGTTCCATAGGTCTTGATTGTCCTGGTGGCCCTGAACCCCAGGTTTTGTCTCCAGTGAGATGCCTGGCCCGGCTTTCTGTGTGACCT";
  
  /*bw.querySequence("test", seq, brv);

  std::cerr << " brv.length " << brv.size() << std::endl;
  for (auto& i : brv) 
    std::cerr << i << std::endl;
  */

  seq = "TGGTGCAACTTAAATCATACTTGTTATATTTGAAAAAGACTGAATATCTCTACAGAAGAAAAAGGATAGAAATAATTCACACCTAACACCTGTAAGTTGTT";
  std::string seq2 = seq + "CGTCTCCACTAAAAATACAAAAATTAAGCCGGGCATGGTGGCACACGCCTGTAATCCCAGC";
  

  for (int i = 0; i < 100000; ++i) {
    
    bw.querySequence("test2", seq2, brv2);
    
    if (i % 100 == 0)
      std::cerr << i << std::endl;

    brv2.clear();

    //for (auto& i : brv2) 
    //  std::cerr<< i << std::endl;

  }

}
