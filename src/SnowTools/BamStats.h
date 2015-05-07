#ifndef SNOWTOOLS_BAMSTATS_H__
#define SNOWTOOLS_BAMSTATS_H__

#include <unordered_map>
#include <cstdint>
#include <iostream>

#include "SnowTools/Histogram.h"
#include "SnowTools/HTSTools.h"

namespace SnowTools{

class BamReadGroup {

  friend class BamStats;
  
 public:

  BamReadGroup() {}

  BamReadGroup(const std::string& name);

  friend ostream& operator<<(std::ostream& out, const BamReadGroup& rg);

  void addRead(Read &r);

 private:

  size_t reads;
  size_t supp;
  size_t unmap;  
  size_t qcfail;
  size_t duplicate;

  Histogram mapq;
  Histogram nm;
  Histogram isize;
  Histogram clip;
  Histogram phred;
  Histogram len;

  std::string m_name;

};

class BamStats
{

 public:
  
  friend std::ostream& operator<<(std::ostream& out, const BamStats& qc);

  void addRead(Read &r);

 private:
  
  std::unordered_map<std::string, BamReadGroup> m_group_map;

};

}

#endif
