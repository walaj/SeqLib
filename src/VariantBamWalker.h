#ifndef VARIANT_VARIANT_BAM_WALKER_H__
#define VARIANT_VARIANT_BAM_WALKER_H__

#include "SnowTools/BamWalker.h"
#include "SnowTools/BamStats.h"

class VariantBamWalker: public SnowTools::BamWalker
{
 public:

  VariantBamWalker() {}

  VariantBamWalker(const std::string in) : SnowTools::BamWalker(in) {}

  void writeVariantBam();
  
  void TrackSeenRead(Read &r);
  
  void printMessage(const SnowTools::ReadCount &rc_main, const Read &r) const;

 private:
  
  SnowTools::BamStats m_stats;

};
#endif
