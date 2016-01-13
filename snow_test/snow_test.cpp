#include <climits>
#include <boost/test/unit_test.hpp>

#include "SnowTools/MiniRules.h"
#include "SnowTools/GenomicRegion.h"
#include "SnowTools/BamWalker.h"

BOOST_AUTO_TEST_CASE( genomic_region_constructors ) {

  // GenomicRegion Constructors
  SnowTools::GenomicRegion gr(0, 0, 10, '+');
  BOOST_CHECK_EQUAL(gr.width(), 11);

  SnowTools::GenomicRegion gr_empty;
  BOOST_TEST(gr_empty.isEmpty());

  SnowTools::GenomicRegion gr2("chrX", "0", "10");
  BOOST_CHECK_EQUAL(gr2.width(), 11);
  BOOST_CHECK_EQUAL(gr2.chr, 22);

  SnowTools::GenomicRegion gr3("X", "0", "10");
  BOOST_TEST(gr2 == gr3);

  BOOST_CHECK_EQUAL(gr.distance(gr2), -1);
  BOOST_CHECK_EQUAL(gr2.distance(gr), -1);

  // check negative inputs
  SnowTools::GenomicRegion grn(-1,-11,-10);
  BOOST_CHECK_EQUAL(grn.chr, -1);
  BOOST_CHECK_EQUAL(grn.pos1, -11);
  BOOST_CHECK_EQUAL(grn.pos2, -10);

  // check strand constructions
  SnowTools::GenomicRegion gra(0,0,0);
  SnowTools::GenomicRegion grb(0,10000,10001, '+');
  SnowTools::GenomicRegion grc(0,0,3, '-');
  BOOST_CHECK_EQUAL(gra.strand, '*');
  BOOST_CHECK_EQUAL(grb.strand, '+');
  BOOST_CHECK_EQUAL(grc.strand, '-');

  // check point string
  BOOST_CHECK_EQUAL(grb.pointString(), "1:10,000(+)");

}

BOOST_AUTO_TEST_CASE( genomic_region_bad_inputs ) {

  BOOST_CHECK_THROW(SnowTools::GenomicRegion(0, 10, 9), std::invalid_argument);

  BOOST_CHECK_THROW(SnowTools::GenomicRegion::chrToString(-1), std::invalid_argument);

  BOOST_CHECK_THROW(SnowTools::GenomicRegion(0,0,0,'P'), std::invalid_argument);

}

BOOST_AUTO_TEST_CASE( genomic_region_random ) {

  SnowTools::GenomicRegion gr; 
  gr.random(42);
  BOOST_CHECK_EQUAL(gr.pointString(), "9:69,477,830(*)");
  
}


BOOST_AUTO_TEST_CASE( genomic_region_check_to_string ) {

  SnowTools::GenomicRegion gr("X", "0","1000");
  BOOST_CHECK_EQUAL(gr.toString(), "X:0-1,000(*)");

  SnowTools::GenomicRegion g2(0, 1, 10, '-');
  BOOST_CHECK_EQUAL(g2.toString(), "1:1-10(-)");
}

BOOST_AUTO_TEST_CASE( genomic_region_constructors_with_headers ) {

  // grab a header 
  SnowTools::BamWalker bw("small.bam");

  // check that it sets the chr number correctly
  SnowTools::GenomicRegion grh("GL000207.1", "0", "10", bw.header());
  BOOST_CHECK_EQUAL(grh.chr, 25);

  // and that it can query the header
  BOOST_CHECK_EQUAL(grh.ChrName(bw.header()), "GL000207.1");

  // check for samtools string
  
}

BOOST_AUTO_TEST_CASE( genomic_check_overlaps ) {

  SnowTools::GenomicRegion gr1(0, 0, 10, '+');
  SnowTools::GenomicRegion gr2(1, 0, 10, '+');

  SnowTools::GenomicRegion gr3(0, 10, 20, '+');
  SnowTools::GenomicRegion gr4(1, 4, 10, '+');

  SnowTools::GenomicRegion gr5(1, 11, 12, '+');

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
