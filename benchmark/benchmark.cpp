//#define USE_BOOST

#ifdef USE_BOOST
#include <boost/timer/timer.hpp>
#endif

#include <cmath>

#define RUN_SEQAN 1
//#define RUN_BAMTOOLS 1
//#define RUN_SEQLIB 1

#ifdef RUN_SEQAN
#include <seqan/bam_io.h>
#endif

#ifdef RUN_SEQLIB
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#endif

#ifdef RUN_BAMTOOLS
#include "api/BamReader.h"
#endif

int main()
{
  
  const size_t limit = 1000000;
  size_t count = 0;

  //std::string bam  = "/xchip/gistic/Jeremiah/GIT/SeqLib/seq_test/test_data/small.bam";
  std::string bam = "/broad/broadsv/NA12878/20120117_ceu_trio_b37_decoy/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam";
  std::string obam = "/xchip/gistic/Jeremiah/GIT/SeqLib/seq_test/tmp_out.bam";

#ifdef USE_BOOST
  boost::timer::auto_cpu_timer t;
#endif

#ifdef RUN_BAMTOOLS
  std::cerr << " **** RUNNING BAMTOOLS **** " << std::endl;
  BamTools::BamReader btr;
  btr.Open(bam);
  BamTools::BamAlignment ba;
  std::vector<BamTools::BamAlignment> bav;
  while(btr.GetNextAlignment(ba) && count++ < limit)
    bav.push_back(ba);
#endif    

#ifdef RUN_SEQLIB
  std::cerr << " **** RUNNING SEQLIB **** " << std::endl;
  SeqLib::BamReader r(bam);
  //SeqLib::BamWriter w(SeqLib::BAM);
  //w.SetHeader(r.Header());
  //w.Open(obam);
  SeqLib::BamRecord rec;
  bool rule = false;
  SeqLib::BamRecordVector brv;
  while(r.GetNextRead(rec, rule) && count++ < limit) 
    brv.push_back(rec);
    //w.writeAlignment(rec);
  return 0;
#endif

#ifdef RUN_SEQAN

  seqan::BamFileIn bamFileIn;
  seqan::BamHeader header;

  if (!open(bamFileIn, bam.c_str(), seqan::OPEN_RDONLY))
    {
      std::cerr << "ERROR: could not open input file " << bam << ".\n";
      return 1;
    }

  /*
  // Open output SAM file.
  seqan::BamFileOut samFileOut(context(bamFileIn), toCString(samFileOutName));
  
  // Copy header.
  try
    {
      readHeader(header, bamFileIn);
      writeHeader(samFileOut, header);
    }
  catch (seqan::IOError const & e)
    {
      std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
    }
  */

  // Copy all records.
  std::vector<seqan::BamAlignmentRecord> record_vec;
  seqan::BamAlignmentRecord record;
  while (!atEnd(bamFileIn) && count++ < limit)
    {
      try
	{
	  readRecord(record, bamFileIn);
	  record_vec.push_back(record);
	  //writeRecord(samFileOut, record);
	}
      catch (seqan::IOError const & e)
	{
	  std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
	}
    }
  
  return 0;

#endif

}
