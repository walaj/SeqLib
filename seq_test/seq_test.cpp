#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include<boost/test/unit_test.hpp>

#include <climits>
#include <boost/test/unit_test.hpp>

#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/ReadFilter.h"
#include "SeqLib/FermiAssembler.h"

#define SBAM "test_data/small.bam"
#define OBAM "test_data/small_out.bam"
#define OCRAM "test_data/small_out.cram"
#define HGREF "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
#define TREF "test_data/test_ref.fa"
#define OREF "tmp_output.fa"

BOOST_AUTO_TEST_CASE( read_filter_0 ) {

  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  // a BAM header check
  BOOST_CHECK_EQUAL(h.GetSequenceLength(0), 249250621);
  BOOST_CHECK_EQUAL(h.GetSequenceLength(3), 191154276);
  BOOST_CHECK_EQUAL(h.GetSequenceLength("1"), 249250621);
  BOOST_CHECK_EQUAL(h.GetSequenceLength("4"), 191154276);
  BOOST_CHECK_EQUAL(h.GetSequenceLength("d4"), -1);
  BOOST_CHECK_EQUAL(h.GetSequenceLength(10000), -1);

  SeqLib::BamRecord rec;
  size_t count = 0;
  
  while(br.GetNextRecord(rec) && count++ < 10000) {
    
  }

}

BOOST_AUTO_TEST_CASE( merge ) {

  SeqLib::GRC grc;
  // add two more that we know of
  grc.add(SeqLib::GenomicRegion(23, 10,100));
  grc.add(SeqLib::GenomicRegion(23, 20,110));

  grc.add(SeqLib::GenomicRegion(2, 10,100));
  grc.add(SeqLib::GenomicRegion(2, 20,110));
  grc.add(SeqLib::GenomicRegion(2, 200,310));

  grc.mergeOverlappingIntervals();
  BOOST_CHECK_EQUAL(grc.size(), 3);
  BOOST_CHECK_EQUAL(grc[0].chr, 2);
  BOOST_CHECK_EQUAL(grc[1].chr, 2);
  BOOST_CHECK_EQUAL(grc[2].chr, 23);
  BOOST_CHECK_EQUAL(grc[2].pos2, 110);
  BOOST_CHECK_EQUAL(grc[2].pos1, 10);
}

BOOST_AUTO_TEST_CASE ( interval_queries ) {

  SeqLib::GRC grc;

  // create a large GRC
  for (int i = 0; i < 10; ++i) {
    int chr = rand() % 23;
    int pos = rand() % 10000;
    grc.add(SeqLib::GenomicRegion(chr, pos, pos + 100));
  }
  grc.mergeOverlappingIntervals();

  // add two more that we know of
  grc.add(SeqLib::GenomicRegion(23, 10,100));
  grc.add(SeqLib::GenomicRegion(23, 20,110));

  // create the interval tree
  grc.createTreeMap();

  SeqLib::OverlapResult orl;
  SeqLib::GRC results = grc.findOverlaps(SeqLib::GenomicRegion(23, 10, 100), true, orl);

  for (auto& i : results)
    std::cerr << " GRC overlaps results " << i << std::endl;
  
  BOOST_CHECK_EQUAL(results.size(), 2);
  BOOST_CHECK_EQUAL(results[1].pos2, 100);

  grc.mergeOverlappingIntervals();

  for(auto& r : grc)
    std::cerr << r << std::endl;

  std::vector<int32_t> q, s;
  results = grc.findOverlaps(grc, q, s, true);

  std::cerr << " results.size " << results.size() << " Input size " << grc.size() << std::endl;
  BOOST_CHECK_EQUAL(results.size(), grc.size());
  BOOST_CHECK_EQUAL(results.width(), grc.width());
  
}

BOOST_AUTO_TEST_CASE( json_parse ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  
  std::string rules = "{\"global\" : {\"!anyflag\" : 1536, \"phred\" : 4}, \"\" : { \"rules\" : [{\"ic\" : true}, {\"clip\" : 5}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";  

  SeqLib::ReadFilterCollection rfc(rules, br.Header());

  std::cerr << rfc << std::endl;

  SeqLib::BamRecord rec;
  size_t count = 0;

  while(br.GetNextRecord(rec) && count++ < 10000) {
    if (!rfc.isValid(rec))
      continue;
    // test global flag rule
    if ( (rec.QCFailFlag() || rec.DuplicateFlag())) {
      std::cerr << rec << std::endl;
      assert(false);
    }
  }

  /// direct from string
  SeqLib::ReadFilterCollection rfc2(rules, br.Header());

  while(br.GetNextRecord(rec) && count++ < 10000) {
    if (!rfc.isValid(rec))
      continue;
    // test global flag rule
    if ( (rec.QCFailFlag() || rec.DuplicateFlag())) {
      std::cerr << rec << std::endl;
      assert(false);
    }
  }

    
}

BOOST_AUTO_TEST_CASE( sw_alignment ) {

  const std::string ref = "ACTGCGAGCGACTAGCTCGTAGCTAGCTAGCTAGCTAGTGACTGCGGGCGATCATCGATCTTTTATTATCGCGATCGCTACGAC";
  const std::string seq = "ACTGCGAGCGACTAGCTCGTAGCTAGCTAGCTAGCTAGTGACTGCGGGCGATCATCGATCTTTTATTATCGCGATCGCTACGAC";
  //const std::string seq =                "CTCGTAGCTAGCTGCTAGCTAGTGACTGCGGGCGATCATCGATCTTTTATTATCGCG";
  const SeqLib::GenomicRegion gr(0,0,0);
  SeqLib::BamRecord b("test_name", seq, ref, &gr);
  
  std::cerr << " SMITH WATERMAN " << std::endl;
  std::cerr << b << std::endl;
}

BOOST_AUTO_TEST_CASE( read_filter_1 ) {

  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  SeqLib::GRC g;
  g.add(SeqLib::GenomicRegion(h.Name2ID("X"), 1100000, 1800000));
  g.createTreeMap();

  // make a new rule set
  SeqLib::ReadFilterCollection rfc;

  // make a new filter region
  SeqLib::ReadFilter rf;

  // add an isize rule on whole-genome
  SeqLib::AbstractRule ar;
  ar.isize = SeqLib::Range(200, 600, false); // 200 to 600, not inverted
  ar.mapq  = SeqLib::Range(10, 50, false); // 200 to 600, not inverted
  ar.nm    = SeqLib::Range(1, 1, false); // 200 to 600, not inverted
  rf.AddRule(ar);

  rf.setRegions(g);

  // add to the filter collection
  rfc.AddReadFilter(rf);

  SeqLib::GRC gback = rfc.getAllRegions();
  BOOST_CHECK_EQUAL(gback.size(), g.size());
  for (size_t i = 0; i < gback.size(); ++i)
    assert(g[i] == gback[i]);
  

  // display
  std::cerr << br.PrintRegions() << std::endl;

  // read / filter the reads
  SeqLib::BamRecord rec;
  size_t count = 0;

  while(br.GetNextRecord(rec) && count++ < 10000) {

    if (!rfc.isValid(rec))
      continue;

    // test isize rule
    if (!(rec.FullInsertSize() >= 200 || rec.FullInsertSize() <= 600)) {
      std::cerr << rec.FullInsertSize() << std::endl;
      assert(false);
    }
    // test mapq rule
    if (!(rec.MapQuality() >= 10 || rec.MapQuality() <= 50)) {
      std::cerr << rec.MapQuality() << std::endl;
      assert(false);
    }
    // test nm rule
    if (!(rec.GetIntTag("NM") != 1)) {
      std::cerr << rec.GetIntTag("NM") << std::endl;
      assert(false);
    }
    
  }
}

BOOST_AUTO_TEST_CASE ( seq_utils ) {
  
  // add commas
  BOOST_CHECK_EQUAL(SeqLib::AddCommas(1),"1");
  BOOST_CHECK_EQUAL(SeqLib::AddCommas(1000000),"1,000,000");
  
  // percent calc
  BOOST_CHECK_EQUAL(SeqLib::percentCalc(10,100), 10);
  BOOST_CHECK_EQUAL(SeqLib::percentCalc(7,8), 87);
  BOOST_CHECK_EQUAL(SeqLib::percentCalc(9,10), 90);
  BOOST_CHECK_EQUAL(SeqLib::percentCalc(2,3), 66);

  // scrub string
  BOOST_CHECK_EQUAL(SeqLib::scrubString("chr1", "chr"), "1");
  BOOST_CHECK_EQUAL(SeqLib::scrubString("chr1", ""), "chr1");
  BOOST_CHECK_EQUAL(SeqLib::scrubString("chr1", "dd"), "chr1");
  BOOST_CHECK_EQUAL(SeqLib::scrubString("chr1", "1"), "chr");

}

BOOST_AUTO_TEST_CASE( bam_record ) {

  // get a record
  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamRecord r;
  
  SeqLib::BamRecordVector brv;
  
  size_t count = 0;
  br.GetNextRecord(r);
  
  BOOST_CHECK_EQUAL(r.asGenomicRegion().chr, 22);
  BOOST_CHECK_EQUAL(r.asGenomicRegion().pos1,999901);
  BOOST_CHECK_EQUAL(r.asGenomicRegion().pos2,1000002);
  BOOST_CHECK_EQUAL(r.asGenomicRegion().strand,'+');

  BOOST_CHECK_EQUAL(r.asGenomicRegionMate().chr, 22);
  BOOST_CHECK_EQUAL(r.asGenomicRegionMate().pos1,999993);
  BOOST_CHECK_EQUAL(r.asGenomicRegionMate().pos2,1000094);
  BOOST_CHECK_EQUAL(r.asGenomicRegionMate().strand,'-');

  BOOST_CHECK_EQUAL(std::floor(r.MeanPhred()), 34);

  BOOST_CHECK_EQUAL(r.CountNBases(), 0);

  r.SetQname("testq");
  BOOST_CHECK_EQUAL(r.Qname(), "testq");

  const std::string s = "ACTGCTAGCTAGCTACTCTGCTACTATATTAGCGCGCATTCGC";
  r.SetSequence(s);
  BOOST_CHECK_EQUAL(r.Sequence(), s);
  
  r.SmartAddTag("ST", "1");
  r.SmartAddTag("ST", "3");
  r.SmartAddTag("ST", "5");
    
  BOOST_CHECK_EQUAL(r.GetSmartIntTag("ST").size(), 3);
  BOOST_CHECK_EQUAL(r.GetSmartIntTag("ST").at(2), 5);
  BOOST_CHECK_EQUAL(r.GetSmartDoubleTag("ST").size(), 3);
  BOOST_CHECK_EQUAL(r.GetSmartDoubleTag("ST").at(2), 5.0);
    
  BOOST_CHECK_EQUAL(r.GetSmartStringTag("ST").size(), 3);
  BOOST_CHECK_EQUAL(r.GetSmartStringTag("ST")[1], "3");
}

BOOST_AUTO_TEST_CASE( fermi_assemble ) {

  SeqLib::FermiAssembler f;

  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamRecord r;

  SeqLib::BamRecordVector brv;
  
  size_t count = 0;
  while(br.GetNextRecord(r) && count++ < 1000) {
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

  std::cout << h.AsString() << std::endl;
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
  
  const std::string in = "2:1,000,000-2,000,000";
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

BOOST_AUTO_TEST_CASE( small_trie_from_file) {

  SeqLib::AbstractRule ar;
  const bool inverted = false;
  ar.addMotifRule("test_data/motif.txt", inverted);

  SeqLib::ReadFilterCollection rfc;
  SeqLib::ReadFilter rf;
  rf.AddRule(ar);
  rfc.AddReadFilter(rf);

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");

  SeqLib::BamRecord rec;
  bool rule;
  size_t count = 0;
  while (br.GetNextRecord(rec) && count++ < 1000){
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

  BOOST_CHECK_THROW(SeqLib::GenomicRegion gr3("X", "a", "10", SeqLib::BamHeader()), std::invalid_argument);
  BOOST_CHECK_THROW(SeqLib::GenomicRegion gr3("X", "1000000000000000000000000000000000000000000000000000000000000000000000000000000", "10", SeqLib::BamHeader()), std::out_of_range);

  BOOST_CHECK_EQUAL(gr.distanceBetweenStarts(gr2), -1);
  BOOST_CHECK_EQUAL(gr2.distanceBetweenStarts(gr), -1);

  SeqLib::BamReader br("test_data/small.bam");
  BOOST_CHECK_EQUAL(SeqLib::GenomicRegion("X","1","100", br.Header()).chr, 22);

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

  // check pretty string
  std::stringstream ss;
  ss << grb;
  BOOST_CHECK_EQUAL(ss.str(), "1:10,000-10,001(+)");

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

  // set some options
  bwa.setGapOpen(32);
  bwa.setGapExtension(1);
  bwa.setMismatchPenalty(18);
  bwa.setAScore(2);
  bwa.setZDropoff(100);
  bwa.set3primeClippingPenalty(5);
  bwa.set5primeClippingPenalty(5);
  bwa.setBandwidth(1000);
  bwa.setReseedTrigger(1.5);

  BOOST_CHECK_THROW(bwa.setGapOpen(-1), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.setGapExtension(-1), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.setMismatchPenalty(-18), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.setAScore(-2), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.setZDropoff(-100), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.set3primeClippingPenalty(-5), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.set5primeClippingPenalty(-5), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.setBandwidth(-1000), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.setReseedTrigger(-1.5), std::invalid_argument);

  // no index loaded yet
  BOOST_CHECK_THROW(bwa.ChrIDToName(1), std::runtime_error);

  // load a test index
  BOOST_TEST(SeqLib::read_access_test(TREF));
  bwa.retrieveIndex(TREF);

  BOOST_CHECK_EQUAL(bwa.NumSequences(), 2);

  BOOST_CHECK_EQUAL(bwa.ChrIDToName(0), "ref1");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(1), "ref2");
  BOOST_CHECK_THROW(bwa.ChrIDToName(2), std::out_of_range);

  BOOST_CHECK(!bwa.retrieveIndex("test_data/small.bam"));

  SeqLib::BamHeader hh = bwa.HeaderFromIndex();
  BOOST_CHECK_EQUAL(hh.NumSequences(), 2);

  // error check the index construction
  SeqLib::USeqVector usv_bad1 = {
    {"ref1", "ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCGATCGATCGATCGTAGC"}, 
    {"ref4", ""},
    {"ref5", "CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT"}};
  SeqLib::USeqVector usv_bad2 = {
    {"", "ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCGATCGATCGATCGTAGC"}, 
    {"ref4", "ACCATCGCAGCAGCTATCTATTATATCGGCAGCATCTAGC"},
    {"ref5", "CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT"}};
  BOOST_CHECK_THROW(bwa.constructIndex(usv_bad1), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.constructIndex(usv_bad2), std::invalid_argument);

  // construct a normal index
  SeqLib::USeqVector usv = {
    {"ref3", "ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCGATCGATCGATCGTAGC"}, 
    {"ref4", "CTACTTTATCATCTACACACTGCCTGACTGCGGCGACGAGCGAGCAGCTACTATCGACT"},
    {"ref5", "CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT"}};


  bwa.constructIndex(usv);

  BOOST_CHECK_EQUAL(bwa.NumSequences(), 3);
  bwa.ChrIDToName(1);

  BOOST_CHECK_THROW(bwa.ChrIDToName(-1), std::out_of_range);
  BOOST_CHECK_THROW(bwa.ChrIDToName(10000), std::out_of_range);

  std::cout << bwa.ChrIDToName(2) << std::endl;

  BOOST_CHECK_EQUAL(bwa.ChrIDToName(0), "ref3");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(1), "ref4");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(2), "ref5");
  BOOST_CHECK_THROW(bwa.ChrIDToName(3), std::out_of_range);


  // write the index
  BOOST_CHECK(bwa.writeIndex(OREF));

  // write the fasta
  std::ofstream os;
  os.open(OREF);
  os << "<" << usv[0].name << std::endl << usv[0].seq <<
    std::endl << usv[1].name << std::endl << usv[1].seq << 
    std::endl << usv[2].name << std::endl << usv[2].seq << 
    std::endl;

  // read it back
  BOOST_CHECK(bwa.retrieveIndex(OREF));

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
  bw.SetRegion(SeqLib::GenomicRegion(22, 1000000, 1001000));

  // make a set of locations
  SeqLib::GRC grc;
  grc.add(SeqLib::GenomicRegion(0, 1, 100));
  grc.add(SeqLib::GenomicRegion(1, 1, 100));

  // set regions
  bw.SetMultipleRegions(grc);

  // write index of new bam
  // should print a warning since no write bam is specified
  //bw.BuildIndex();

  // open an output BAM
  //bw.OpenWriteBam(OBAM);

  // set tags to strip
  //bw.setStripTags("OQ,BI");

  // loop through and grab some reads
  SeqLib::BamRecord r;
  size_t count = 0;
  while (bw.GetNextRecord(r)) {
    //if (++count % 10 == 0)
    //  bw.WriteRecord(r);
  }
  
  // display info about BAM
  std::cerr << bw << std::endl;

  // write index of new bam
  //bw.BuildIndex();

  // reset the walker
  bw.Reset();

  // write as a cram
  //bw.OpenWriteBam(OCRAM);
    
  //
  //bw.setCram(OCRAM, HGREF);

  // print cram writer
  //std::cerr << bw << std::endl;
  // write the CRAM
  //while (bw.GetNextRecord(r, rule)) {
  //  if (++count % 10 == 0) {
  //    std::cerr << count << std::endl;
  //    bw.WriteRecord(r);
  //  }
  //}

}


BOOST_AUTO_TEST_CASE( sequtils ) {

  std::string seq = "actgACGTnTCN";

  SeqLib::rcomplement(seq);
  
  BOOST_CHECK_EQUAL(seq, "NGAnACGTcagt");

}

BOOST_AUTO_TEST_CASE( gr_random ) {

  SeqLib::GenomicRegion gr;
  gr.random();
  std::cerr << " RANDOM " << gr << std::endl;

}

BOOST_AUTO_TEST_CASE( bam_write ) {


  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  SeqLib::BamRecord rec;

  // empty constructor
  SeqLib::BamWriter w;
  
  // specify bam explicitly
  //w = SeqLib::BamWriter(SeqLib::BAM);

  BOOST_CHECK_THROW(w.WriteHeader(), std::runtime_error);
  BOOST_CHECK_THROW(w.Close(), std::runtime_error);
  BOOST_CHECK_THROW(w.BuildIndex(), std::runtime_error);
  BOOST_CHECK_THROW(w.WriteRecord(rec), std::runtime_error);

  w.Open("tmp_out.bam");

  BOOST_CHECK_THROW(w.WriteHeader(), std::runtime_error);

  w.SetHeader(h);

  w.WriteHeader();

  size_t count = 0;

  while(br.GetNextRecord(rec) && count++ < 10000) {
    w.WriteRecord(rec);
  }


  BOOST_CHECK_THROW(w.BuildIndex(), std::runtime_error);
  w.Close();

  w.BuildIndex();

  // print some info
  std::cerr << w << std::endl;

}

BOOST_AUTO_TEST_CASE( bam_record_more ) {

  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  SeqLib::BamRecord rec;
  size_t count = 0;
  
  while(br.GetNextRecord(rec) && count++ < 100) {
    rec.ClearSeqQualAndTags();
    assert(rec.Sequence().empty());
    assert(rec.Qualities().empty());
    assert(rec.GetIntTag("NM") == 0);
    assert(rec.GetZTag("XA").empty() == (rec.CountBWASecondaryAlignments()==0));
    rec.CountBWAChimericAlignments();
  }

  br.Reset();

  SeqLib::ReadFilterCollection rf;

}

BOOST_AUTO_TEST_CASE( bam_record_manipulation ) {

  SeqLib::Cigar cig;

  // manually construct a cigar
  cig.add(SeqLib::CigarField('M', 10));
  cig.add(SeqLib::CigarField('I', 1));
  cig.add(SeqLib::CigarField('M', 10));
  cig.add(SeqLib::CigarField('D', 1));
  cig.add(SeqLib::CigarField('M', 10));
  cig.add(SeqLib::CigarField('S', 10));
  
  // make a sequence
  const std::string seq = std::string(10, 'A') + std::string(1, 'T') + std::string(10, 'C') + std::string(10, 'G') + std::string(10, 'A');

  // check   
  BOOST_CHECK_EQUAL(cig.NumQueryConsumed(), 41);
  BOOST_CHECK_EQUAL(cig.NumReferenceConsumed(), 31);

  std::stringstream ss;
  ss << cig;

  // cigar from string
  SeqLib::Cigar cig2 = SeqLib::cigarFromString(ss.str());

  // check that the string from / to are consistent
  assert(cig == cig2);
  assert(!(cig != cig2));
  for (int i = 0; i < cig.size(); ++i)
    assert(cig[i] == cig2[i]);
  for (int i = 0; i < cig.size(); ++i)
    assert(!(cig[i] != cig2[i]));

  // manually make a read
  SeqLib::GenomicRegion gr_wrong(0, 100, 131); 
  SeqLib::GenomicRegion gr(0, 100, 130); 
  
  BOOST_CHECK_THROW(SeqLib::BamRecord("dumname", seq, &gr_wrong, cig), std::invalid_argument);
  BOOST_CHECK_THROW(SeqLib::BamRecord("dumname", seq + "A", &gr, cig), std::invalid_argument);

  SeqLib::BamRecord br("dumname", seq, &gr, cig);

  BOOST_CHECK_EQUAL(br.Sequence(), seq);
  BOOST_CHECK_EQUAL(br.GetCigar(), cig);
  BOOST_CHECK_EQUAL(br.Qname(), "dumname");
  BOOST_CHECK_EQUAL(br.Position(), 100);
  BOOST_CHECK_EQUAL(br.Length(), 41);
  BOOST_CHECK_EQUAL(br.ChrID(), 0);

}

BOOST_AUTO_TEST_CASE( change_bam_record ) {

  // get a record
  SeqLib::BamReader br("test_data/small.bam");
  SeqLib::BamRecord r;
  
  SeqLib::BamRecordVector brv;
  
  size_t count = 0;
  br.GetNextRecord(r);

  SeqLib::Cigar c = r.GetCigar();
  std::cerr << c << std::endl;

  // try replace with cigar of same size
  SeqLib::Cigar c2;
  c2.add(SeqLib::CigarField('S', 101));
  r.SetCigar(c2);
  std::cerr << r << std::endl;

  // try replace with new cigar
  SeqLib::Cigar c3;
  c3.add(SeqLib::CigarField('S', 10));
  c3.add(SeqLib::CigarField('M', 91));
  r.SetCigar(c3);
  std::cerr << r << std::endl;

  const std::string new_seq = "ACTGGACTACAC";

  r.SetSequence(new_seq);
  std::cerr << r << std::endl;

  r.SetQname("dummy_qname");
  std::cerr << r << std::endl;

}
