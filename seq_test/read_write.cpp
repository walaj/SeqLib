#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include<boost/test/unit_test.hpp>

#include <climits>
#include <boost/test/unit_test.hpp>

#include "SeqLib/BamWriter.h"
#include "SeqLib/BamReader.h"

BOOST_AUTO_TEST_CASE( stdinput ) {

  /*
  // read a BAM from stdin
  SeqLib::BamReader b("-"); 

  // write it back out
  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    std::cerr << " STDIN " << r << std::endl;
  }
  */
}

BOOST_AUTO_TEST_CASE( cramin ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.cram"); 

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    std::cerr << "CRAM " << r << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( cramin_new_ref ) {

  // read a BAM from stdin
  SeqLib::BamReader b;
  b.OpenReadBam("test_data/small.cram");

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    std::cerr << "CRAM " << r << std::endl;
  }
}


BOOST_AUTO_TEST_CASE( bamin ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.bam"); 

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    std::cerr << "BAM " << r << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( samin ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    std::cerr << "SAM " << r << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( bamout ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::BAM);
  //SeqLib::BamWriter w;
  w.Open("tmp_out.bam");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    w.writeAlignment(r);
  }
  w.CloseBam();
  
}

BOOST_AUTO_TEST_CASE( samout ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::SAM);
  w.Open("tmp_out.sam");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    w.writeAlignment(r);
  }
  w.CloseBam();
  
}


BOOST_AUTO_TEST_CASE( cramout ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.cram"); 

  SeqLib::BamWriter w(SeqLib::CRAM);
  w.Open("tmp_out.cram");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    w.writeAlignment(r);
  }
  w.CloseBam();
  
}

BOOST_AUTO_TEST_CASE( samout_to_stdout ) {

  /*
  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::SAM);
  w.Open("-");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    w.writeAlignment(r);
  }
  w.CloseBam();
  */
}

BOOST_AUTO_TEST_CASE( bamout_to_stdout ) {

  /**
   * dont actually run every time
   * too much stdout-ing
   */

  /*
  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::BAM);
  w.Open("-");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  bool rule;
  size_t count = 0;
  while(b.GetNextRead(r, rule) && count++ < 1) {
    w.writeAlignment(r);
  }
  w.CloseBam();

  */
  
}
