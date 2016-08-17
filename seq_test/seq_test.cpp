#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include<boost/test/unit_test.hpp>

#include <climits>
#include <boost/test/unit_test.hpp>

#include "SeqLib/GenomicRegion.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BamWalker.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/ReadFilter.h"
#include "SeqLib/FermiAssembler.h"

#define SBAM "test_data/small.bam"
#define OBAM "test_data/small_out.bam"
#define OCRAM "test_data/small_out.cram"
#define HGREF "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
#define TREF "test_data/test_ref.fa"
#define OREF "tmp_output.fa"

BOOST_AUTO_TEST_CASE( fermi_assemble ) {

  SeqLib::FermiAssembler f;

  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamRecord r;
  bool rule;
  
  SeqLib::BamRecordVector brv;
  
  size_t count = 0;
  while(br.GetNextRead(r, rule) && count++ < 1000) {
    brv.push_back(r);
  }
  
  f.AddReads(brv);

  f.CorrectReads();

  SeqLib::UnalignedReadVector reads = f.GetSequences();
  BOOST_CHECK_EQUAL(reads.size(), brv.size());

  for (int i = 0; i < reads.size(); ++i) {
    if (brv[i].Sequence() != reads[i].seq) {
      std::cerr << "************" << std::endl;
      std::cerr << brv[i].Sequence() << std::endl;
      std::cerr << reads[i].seq << std::endl;
    }
  }
    
  // peform the assembly
  std::cerr << "...performing assembly" << std::endl;
  f.PerformAssembly();

  // retrieve the contigs
  std::vector<std::string> contigs = f.GetContigs();

}

BOOST_AUTO_TEST_CASE( bam_header_stdout ) {

  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();
  
  h.WriteToStdout();
  
}

BOOST_AUTO_TEST_CASE( bam_header_name2id ) {

  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  BOOST_CHECK_EQUAL(h.Name2ID("2"), 1);  
  BOOST_CHECK_EQUAL(h.Name2ID("23"), -1);  

}

BOOST_AUTO_TEST_CASE( bam_header_id2name ) {
  
  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();
  
  BOOST_CHECK_EQUAL(h.IDtoName(2), "3");
  BOOST_CHECK_THROW(h.IDtoName(100), std::out_of_range);
  BOOST_CHECK_THROW(h.IDtoName(-1), std::invalid_argument);
  BOOST_CHECK_THROW(SeqLib::BamHeader().IDtoName(1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE( genomic_ranges_string_constructor) {
  
  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();
  
  const std::string in = "chr2:1,000,000-2,000,000";
  SeqLib::GenomicRegion gr(in, h);
  BOOST_CHECK_EQUAL(gr.chr, 1);
  BOOST_CHECK_EQUAL(gr.pos1, 1000000);
  BOOST_CHECK_EQUAL(gr.pos2, 2000000);

  BOOST_CHECK_THROW(SeqLib::GenomicRegion(in, SeqLib::BamHeader()), std::invalid_argument);

  BOOST_CHECK_EQUAL(gr.ChrName(h), "2");
  BOOST_CHECK_EQUAL(gr.ChrName(SeqLib::BamHeader()), "2");
  gr.chr = 1000;
  BOOST_CHECK_THROW(gr.ChrName(h), std::invalid_argument);

}

BOOST_AUTO_TEST_CASE( genomic_region_less_than ) {

  SeqLib::GenomicRegion gr1(0, 1, 2);
  SeqLib::GenomicRegion gr2(1, 1, 2);
  SeqLib::GenomicRegion gr3(1, 2, 2);
  SeqLib::GenomicRegion gr4(1, 6, 6);

  BOOST_CHECK(gr1 < gr2);
  BOOST_CHECK(gr2 > gr1);
  BOOST_CHECK(!(gr1 > gr2));

  BOOST_CHECK(gr2 < gr3);
  BOOST_CHECK(gr3 > gr2);
  BOOST_CHECK(!(gr2 > gr3));

  BOOST_CHECK(gr3 < gr4);
  BOOST_CHECK(!(gr4 == gr3));
  BOOST_CHECK(!(gr3 > gr4));
  BOOST_CHECK(gr4 > gr3);

  BOOST_CHECK(!(gr1 < gr1));
  BOOST_CHECK(!(gr1 > gr1));

  BOOST_CHECK(!(gr1 != gr1));
  BOOST_CHECK(gr2 != gr1);
  BOOST_CHECK(gr3 != gr1);
  BOOST_CHECK(gr4 != gr3);

  BOOST_CHECK(gr1 >= gr1);
  BOOST_CHECK(gr2 >= gr2);
  BOOST_CHECK(gr3 >= gr3);
  BOOST_CHECK(gr4 >= gr4);

  BOOST_CHECK(gr1 <= gr1);
  BOOST_CHECK(gr2 <= gr2);
  BOOST_CHECK(gr3 <= gr3);
  BOOST_CHECK(gr4 <= gr4);

  BOOST_CHECK(gr1 <= gr2);
  BOOST_CHECK(gr2 >= gr1);

  BOOST_CHECK(gr2 <= gr3);
  BOOST_CHECK(gr3 >= gr2);

}

BOOST_AUTO_TEST_CASE( genomic_region_distance ) {

  SeqLib::GenomicRegion gr1(0, 10, 100);
  SeqLib::GenomicRegion gr2(0, 10, 200);
  SeqLib::GenomicRegion gr3(1, 10, 100);
  SeqLib::GenomicRegion gr4(0, 100, 100);

  BOOST_CHECK_EQUAL(gr1.distanceBetweenEnds(gr3), -1);
  BOOST_CHECK_EQUAL(gr1.distanceBetweenEnds(gr1), 0);
  BOOST_CHECK_EQUAL(gr1.distanceBetweenEnds(gr2), 100);
  BOOST_CHECK_EQUAL(gr1.distanceBetweenEnds(gr4), 0);

  BOOST_CHECK_EQUAL(gr1.distanceBetweenStarts(gr3), -1);
  BOOST_CHECK_EQUAL(gr1.distanceBetweenStarts(gr1), 0);
  BOOST_CHECK_EQUAL(gr1.distanceBetweenStarts(gr2), 0);
  BOOST_CHECK_EQUAL(gr1.distanceBetweenStarts(gr4), 90);

}

BOOST_AUTO_TEST_CASE( stdinput ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.bam");

  // write it back out
  SeqLib::BamWriter w(SeqLib::SAM);
  
  w.SetHeader(b.Header());
  w.Open("tmp_out_from_stdin.bam");
  w.WriteHeader();
  
  SeqLib::BamRecord r;
  bool rule;
  while(b.GetNextRead(r, rule)) {
    w.writeAlignment(r);
  }

}

BOOST_AUTO_TEST_CASE( large_trie ) {

  const std::string dictionary = "ACTG";
  
  const int string_size = 20;
  const int string_count = 10000;

  SeqLib::AhoCorasick aho;

  std::vector<std::string> k;

  std::cerr << "...generating key" << std::endl;
  for (int i = 0; i < string_count; ++i) {
    char* c = (char*) malloc(string_size + 1);
    for (int j = 0; j < string_size; ++j)
      c[j] = dictionary.at(rand() % 4);
    c[string_size] = '\0';
    k.push_back(std::string(c));
    free(c);
  }
  std::cerr << "...done with key" << std::endl;

  std::cerr << "...generating trie" << std::endl;
  for (auto& i : k)
    aho.AddMotif(i);
  std::cerr << "...done generating trie" << std::endl;

  std::cerr << "...querying trie" << std::endl;
  auto result = aho.aho_trie->parse_text(k[0]);
  std::cerr << "...querying trie fast" << std::endl;  
  for (int i = 0; i < string_count; ++i) {
    //if (i % 20000 == 0)
    //  std::cerr << "... " << i << std::endl;
    auto result = aho.aho_trie->parse_text(k[i]);
  }
    
}

BOOST_AUTO_TEST_CASE( genomic_region_constructors ) {

  // GenomicRegion Constructors
  SeqLib::GenomicRegion gr(0, 0, 10, '+');
  BOOST_CHECK_EQUAL(gr.width(), 11);

  SeqLib::GenomicRegion gr_empty;
  BOOST_TEST(gr_empty.isEmpty());

  SeqLib::GenomicRegion gr2("chrX", "0", "10", SeqLib::BamHeader());
  BOOST_CHECK_EQUAL(gr2.width(), 11);
  BOOST_CHECK_EQUAL(gr2.chr, 22);

  SeqLib::GenomicRegion gr3("X", "0", "10", SeqLib::BamHeader());
  BOOST_TEST(gr2 == gr3);

  BOOST_CHECK_EQUAL(gr.distanceBetweenStarts(gr2), -1);
  BOOST_CHECK_EQUAL(gr2.distanceBetweenStarts(gr), -1);

  // check negative inputs
  SeqLib::GenomicRegion grn(-1,-11,-10);
  BOOST_CHECK_EQUAL(grn.chr, -1);
  BOOST_CHECK_EQUAL(grn.pos1, -11);
  BOOST_CHECK_EQUAL(grn.pos2, -10);

  // check strand constructions
  SeqLib::GenomicRegion gra(0,0,0);
  SeqLib::GenomicRegion grb(0,10000,10001, '+');
  SeqLib::GenomicRegion grc(0,0,3, '-');
  BOOST_CHECK_EQUAL(gra.strand, '*');
  BOOST_CHECK_EQUAL(grb.strand, '+');
  BOOST_CHECK_EQUAL(grc.strand, '-');

  // check point string
  BOOST_CHECK_EQUAL(grb.pointString(), "1:10,000(+)");

}

BOOST_AUTO_TEST_CASE( genomic_region_bad_inputs ) {

  BOOST_CHECK_THROW(SeqLib::GenomicRegion(0, 10, 9), std::invalid_argument);

  //BOOST_CHECK_THROW(SeqLib::GenomicRegion::chrToString(-1), std::invalid_argument);

  BOOST_CHECK_THROW(SeqLib::GenomicRegion(0,0,0,'P'), std::invalid_argument);

}

// BOOST_AUTO_TEST_CASE( genomic_region_random ) {

//   SeqLib::GenomicRegion gr; 
//   std::srand(42);
//   gr.random();
//   BOOST_CHECK_EQUAL(gr.pointString(), "9:69,477,830(*)");
  
// }

BOOST_AUTO_TEST_CASE( genomic_region_range_operations ) {

  SeqLib::GenomicRegion gr(0,1,10);
  SeqLib::GenomicRegion gr2(0,1,11);
  gr.pad(3);
  gr2.pad(-3);
  BOOST_CHECK_EQUAL(gr.pos1,-2);
  BOOST_CHECK_EQUAL(gr.pos2,13);
  BOOST_CHECK_EQUAL(gr2.pos1,4);
  BOOST_CHECK_EQUAL(gr2.pos2,8);

  BOOST_CHECK_THROW(gr.pad(-10), std::out_of_range);

}

BOOST_AUTO_TEST_CASE( genomic_region_check_to_string ) {

  SeqLib::GenomicRegion gr("X", "0","1000", SeqLib::BamHeader());
  BOOST_CHECK_EQUAL(gr.toString(), "X:0-1,000(*)");

  SeqLib::GenomicRegion g2(0, 1, 10, '-');
  BOOST_CHECK_EQUAL(g2.toString(), "1:1-10(-)");

  // check default ref to string conversion (no header)
  BOOST_CHECK_EQUAL(gr.ChrName(SeqLib::BamHeader()), "X");
}

BOOST_AUTO_TEST_CASE( genomic_check_overlaps ) {

  SeqLib::GenomicRegion gr1(0, 0, 10, '+');
  SeqLib::GenomicRegion gr2(1, 0, 10, '+');

  SeqLib::GenomicRegion gr3(0, 10, 20, '+');
  SeqLib::GenomicRegion gr4(1, 4, 10, '+');

  SeqLib::GenomicRegion gr5(1, 11, 12, '+');

  // partial overlaps should be one
  BOOST_CHECK_EQUAL(gr1.getOverlap(gr3), 1);

  // argument contained gets 2
  BOOST_CHECK_EQUAL(gr2.getOverlap(gr4), 2);

  // object contained gets 3 
  BOOST_CHECK_EQUAL(gr4.getOverlap(gr2), 3);

  // same chr, no overlap
  BOOST_CHECK_EQUAL(gr4.getOverlap(gr5), 0);
  BOOST_CHECK_EQUAL(gr5.getOverlap(gr4), 0);

}

BOOST_AUTO_TEST_CASE( bwa_wrapper ) {

  SeqLib::BWAWrapper bwa;

  // load a test index
  BOOST_TEST(SeqLib::read_access_test(TREF));
  bwa.retrieveIndex(TREF);

  BOOST_CHECK_EQUAL(bwa.NumSequences(), 2);

  BOOST_CHECK_EQUAL(bwa.ChrIDToName(0), "ref1");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(1), "ref2");
  BOOST_CHECK_THROW(bwa.ChrIDToName(2), std::out_of_range);

  SeqLib::BamHeader hh = bwa.HeaderFromIndex();
  BOOST_CHECK_EQUAL(hh.NumSequences(), 2);

  SeqLib::USeqVector usv = {
    {"ref3", "ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCGATCGATCGATCGTAGC"}, 
    {"ref4", "CTACTTTATCATCTACACACTGCTACTGACTGCGGCGACGAGCGAGCAGCTACTATCGACT"},
    {"ref5", "CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT"}};

  bwa.constructIndex(usv);

  BOOST_CHECK_EQUAL(bwa.ChrIDToName(0), "ref3");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(1), "ref4");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(2), "ref5");
  BOOST_CHECK_THROW(bwa.ChrIDToName(3), std::out_of_range);

  // write the index
  bwa.writeIndex(OREF);

  // write the fasta
  std::ofstream os;
  os.open(OREF);
  os << "<" << usv[0].name << std::endl << usv[0].seq <<
    std::endl << usv[1].name << std::endl << usv[1].seq << 
    std::endl << usv[2].name << std::endl << usv[2].seq << 
    std::endl;

  // read it back
  bwa.retrieveIndex(OREF);

  // check that its good
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(0), "ref3");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(1), "ref4");
  
  // try aligning a sequence
  std::cerr << "...aligning sequences" << std::endl;
  SeqLib::BamRecordVector brv, brv2;
  bool hardclip = false;
  bwa.alignSingleSequence("ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCG", "name", brv, 0.9, hardclip, 1);
  // reverse complement
  bwa.alignSingleSequence("CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGC", "name", brv2, 0.9, hardclip, 2);

  std::cerr << "...checking aligned sequences" << std::endl;
  std::cerr << "...brv.size() " << brv.size() << std::endl;
  BOOST_CHECK_EQUAL(brv[0].Qname(), "name");
  BOOST_CHECK_EQUAL(brv[0].ChrID(), 2);
  BOOST_CHECK_EQUAL(brv[0].Sequence(), "CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT");
  std::cerr << " brv[0].GetCigar() " << brv[0].GetCigar() << std::endl;
  BOOST_CHECK_EQUAL(brv[0].GetCigar()[0].Type(), 'M');
  BOOST_CHECK_EQUAL(brv[0].GetCigar()[0].Length(), 38);

  // check that it got both alignments
  BOOST_CHECK_EQUAL(brv2.size(), 2);

  // print info 
  std::cerr << bwa << std::endl;
  
}

BOOST_AUTO_TEST_CASE( bam_reader ) {

  SeqLib::BamReader bw(SBAM);

  // open index
  bw.setBamReaderRegion(SeqLib::GenomicRegion(22, 1000000, 1001000));

  // make a set of locations
  SeqLib::GRC grc;
  grc.add(SeqLib::GenomicRegion(0, 1, 100));
  grc.add(SeqLib::GenomicRegion(1, 1, 100));

  // set regions
  bw.setBamReaderRegions(grc.asGenomicRegionVector());

  // write index of new bam
  // should print a warning since no write bam is specified
  //bw.makeIndex();

  // open an output BAM
  //bw.OpenWriteBam(OBAM);

  // set tags to strip
  //bw.setStripTags("OQ,BI");

  // loop through and grab some reads
  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while (bw.GetNextRead(r, rule)) {
    //if (++count % 10 == 0)
    //  bw.writeAlignment(r);
  }
  
  // display info about BAM
  std::cerr << bw << std::endl;

  // write index of new bam
  //bw.makeIndex();

  // reset the walker
  bw.resetAll();

  // write as a cram
  //bw.OpenWriteBam(OCRAM);
    
  //
  //bw.setCram(OCRAM, HGREF);

  // print cram writer
  //std::cerr << bw << std::endl;
  // write the CRAM
  //while (bw.GetNextRead(r, rule)) {
  //  if (++count % 10 == 0) {
  //    std::cerr << count << std::endl;
  //    bw.writeAlignment(r);
  //  }
  //}

}


BOOST_AUTO_TEST_CASE( sequtils ) {

  std::string seq = "actgACGTnTCN";

  SeqLib::rcomplement(seq);
  
  BOOST_CHECK_EQUAL(seq, "NGAnACGTcagt");


}
