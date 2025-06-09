#include <catch2/catch.hpp>
#include "SeqLib/BamReader.h"
#include <string>

// convenience macro for the test data path
static const std::string bam_path = std::string(TEST_DATA_DIR) + "/sim.sorted.bam";

TEST_CASE("BamReader opens BAM and reads header", "[BamReader]") {
    SeqLib::BamReader reader;
    REQUIRE(reader.Open(bam_path));            // can open file
    const SeqLib::BamHeader &hdr = reader.Header();
    // our tiny_ref.fa was indexed as a single contig named "BCRABL"
    REQUIRE(hdr.NumSequences() == 4);         
    REQUIRE(hdr.IDtoName(0) == "bcr"); // name matches
    REQUIRE(hdr.IDtoName(1) == "abl"); // name matches    
    reader.Close();
}

TEST_CASE("BamReader iterates through all records", "[BamReader]") {
    SeqLib::BamReader reader;
    REQUIRE(reader.Open(bam_path));

    bool has_supp = false;
    bool has_paired = false;
    bool has_inter = false;
    
    size_t count = 0;
    while (auto rec = reader.Next()) {
        // sanity check: each record should map to one of the references
        REQUIRE(rec->ChrID() < 4);
        ++count;

	if (rec->SupplementaryFlag())
	  has_supp = true;
	if (rec->PairedFlag())
	  has_paired = true;
	if (rec->Interchromosomal())
	  has_inter = true;
    }
    // we expect >0 reads
    REQUIRE(count > 0);

    // expect some supplementary alignments
    REQUIRE(has_supp);
    REQUIRE(has_paired);
    REQUIRE(has_inter);
    
    // once we hit the end, Next() remains nullopt
    REQUIRE_FALSE(reader.Next().has_value());
    reader.Close();
}

TEST_CASE("BamReader Reset allows multiple passes", "[BamReader]") {
    SeqLib::BamReader reader;
    REQUIRE(reader.Open(bam_path));

    // first pass
    size_t pass1 = 0;
    while (reader.Next()) ++pass1;
    REQUIRE(pass1 > 0);

    // reset and do it again
    reader.Reset();
    size_t pass2 = 0;
    while (reader.Next()) ++pass2;
    REQUIRE(pass2 == pass1);

    reader.Close();
}

TEST_CASE("BamReader can seek on regions") {

  SeqLib::BamReader reader;
  REQUIRE(reader.Open(bam_path));

  SeqLib::GRC grc;
  grc.add(SeqLib::GenomicRegion(0,1,1000));
  grc.add(SeqLib::GenomicRegion(1,1,1000));
  grc.add(SeqLib::GenomicRegion(2,1,1000));
  
  reader.SetRegions(grc);

  size_t pass1 = 0;
  while (reader.Next()) ++pass1;
  REQUIRE(pass1 > 1000);

  reader.Close();
}
