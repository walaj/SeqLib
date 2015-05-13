#include "VariantBamWalker.h"

using SnowTools::ReadCount;

void VariantBamWalker::writeVariantBam() 
{

  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);

  ReadCount rc_main;
  Read r;
  std::cout << "...starting run through BAM" << std::endl;
  std::string rule;

  while (GetNextRead(r, rule))
    {
      
      TrackSeenRead(r);
      // read is valid
      if (rule.length()) {
	++rc_main.keep;
	WriteAlignment(r);
      } 

      if (++rc_main.total % 1000000 == 0)
	printMessage(rc_main, r);

    }
}

void VariantBamWalker::TrackSeenRead(Read &r)
{
  m_stats.addRead(r);
}

void VariantBamWalker::printMessage(const ReadCount &rc_main, const Read &r) const 
{
  
  char buffer[100];
  std::string posstring = SnowTools::AddCommas<int>(r_pos(r));
  std::sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
		rc_main.totalString().c_str(), SnowTools::GenomicRegion::chrToString(r_id(r)).c_str(), posstring.c_str(),  
		rc_main.keepString().c_str(), rc_main.percent());
  
  std::printf ("%s | ",buffer);
  SnowTools::displayRuntime(start);
  std::cout << std::endl;
  
}
