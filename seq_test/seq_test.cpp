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
#include "SeqLib/SeqPlot.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BLATWrapper.h"

#define SBAM "test_data/small.bam"
#define OBAM "test_data/small_out.bam"
#define OCRAM "test_data/small_out.cram"
#define HGREF "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
#define TREF "test_data/test_ref.fa"
#define OREF "tmp_output.fa"
#define BEDFILE "test_data/test.bed"
#define VCFFILE "test_data/test.vcf"
#define JSON1 "test_data/example4.json"

using namespace SeqLib;

BOOST_AUTO_TEST_CASE( BLAT_wrapper ) {

  BLATWrapper b;
  b.loadIndex("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta", "/xchip/gistic/Jeremiah/blat/11.ooc");

  SeqLib::BamReader r;
  r.Open(SBAM);
  b.addHeader(r.Header().get_());
  BamRecord rec;
  BamWriter w;
  w.Open("tmp.blat.bam");
  w.SetHeader(r.Header());
  w.WriteHeader();
  size_t count = 0;
  while (r.GetNextRecord(rec) && ++count < 10000) {
    BamRecordVector brv;
    b.querySequence(rec.Qname(), rec.Sequence(), brv);
    for (auto& i : brv)
      w.WriteRecord(i);
  }
  w.Close();

}
