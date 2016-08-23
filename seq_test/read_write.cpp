#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include<boost/test/unit_test.hpp>

#include <climits>
#include <boost/test/unit_test.hpp>

#include "SeqLib/BamWriter.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamPolyReader.h"
#include "SeqLib/SeqPlot.h"

BOOST_AUTO_TEST_CASE( stdinput ) {

#ifdef RUN_STDIN
  // read a BAM from stdin
  SeqLib::BamReader b("-"); 

  // write it back out
  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    std::cerr << " STDIN " << r << std::endl;
  }
#endif
}

BOOST_AUTO_TEST_CASE( cramin ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.cram"); 

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    std::cerr << "CRAM " << r << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( cramin_new_ref ) {

  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.cram");

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 10) {
    std::cerr << "CRAM " << r << std::endl;
  }
}


BOOST_AUTO_TEST_CASE( bamin ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.bam"); 

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    std::cerr << "BAM " << r << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( samin ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
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
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
  
}

BOOST_AUTO_TEST_CASE( samout ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::SAM);
  w.Open("tmp_out.sam");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
  
}


BOOST_AUTO_TEST_CASE( cramout ) {

  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.cram"); 

  SeqLib::BamWriter w(SeqLib::CRAM);
  w.Open("tmp_out.cram");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
  
}

BOOST_AUTO_TEST_CASE( samout_to_stdout ) {

#ifdef RUN_SAM_STDOUT
  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::SAM);
  w.Open("-");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
#endif
}

BOOST_AUTO_TEST_CASE( bamout_to_stdout ) {

  //
  // dont actually run every time
  // too much stdout-ing
  //

#ifdef RUN_BAM_STDOUT
  // read a BAM from stdin
  SeqLib::BamReader b("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::BAM);
  w.Open("-");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
#endif
  
}

BOOST_AUTO_TEST_CASE( bam_poly ) {

  SeqLib::BamPolyReader r;
  
  r.Open("test_data/small.bam");
  r.Open("test_data/small.cram");

  r.SetRegion(SeqLib::GenomicRegion(r.Header().Name2ID("X"),1001000, 1001100));

  SeqLib::BamWriter w(SeqLib::BAM);
  w.Open("tmp_out_poly.bam");
  w.SetHeader(r.Header());
  w.WriteHeader();

  SeqLib::BamRecord rec;
  while(r.GetNextRecord(rec)) {
    w.WriteRecord(rec);
  }

}


BOOST_AUTO_TEST_CASE( plot_test ) {

  SeqLib::BamReader r;
  r.Open("test_data/small.bam");

  SeqLib::GenomicRegion gr("X:1,002,942-1,003,294", r.Header());
  r.SetRegion(gr);

  SeqLib::SeqPlot s;

  s.SetView(gr);

  SeqLib::BamRecord rec;
  SeqLib::BamRecordVector brv;
  while(r.GetNextRecord(rec))
    if (!rec.CountNBases() && rec.MappedFlag())
      brv.push_back(rec);

  s.SetPadding(20);

  std::cout << s.PlotAlignmentRecords(brv);

}
