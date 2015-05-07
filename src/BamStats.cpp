#include "SnowTools/BamStats.h"

#include <cmath>

//#define DEBUG_STATS 1
namespace SnowTools {

BamReadGroup::BamReadGroup(const std::string& name) : reads(0), supp(0), unmap(0), qcfail(0), 
						      duplicate(0), m_name(name)
{

  mapq = Histogram(0,100,1);
  nm = Histogram(0,100,1);
  isize = Histogram(-2,2000,10);
  clip = Histogram(0,100,5);
  phred = Histogram(0,100,1);
  len = Histogram(0,250,1);

}

void BamReadGroup::addRead(Read &r)
{
  int mapqr = r_mapq(r);
  if (mapqr >=0 && mapqr <= 100)
    mapq.addElem(mapqr);
  
  int32_t this_nm;
  r_get_int32_tag(r, "NM", this_nm);
  if (this_nm <= 100)
    nm.addElem(this_nm);
  
  int32_t isizer = -1;
  if (!r_is_pmapped(r))
    isizer = -2;
  else if (r_id(r) == r_mid(r))
    isizer = std::abs(r_isize(r));
  isize.addElem(isizer);

  int32_t c;
  r_get_clip(r,c);
  clip.addElem(c);
  
  len.addElem(r_length(r));

}

void BamStats::addRead(Read &r)
{

  // get the read group
  std::string rg;
  r_get_Z_tag(r, "RG", rg);

#ifdef DEBUG_STATS
  std::cout << "got read group tag " << rg << std::endl;
#endif
  
  std::unordered_map<std::string, BamReadGroup>::iterator ff = m_group_map.find(rg);
  
  if (ff == m_group_map.end())
    {
      m_group_map[rg] = BamReadGroup(rg);
      m_group_map[rg].addRead(r);
    } 
  else
    {
      ff->second.addRead(r);
    }
}

}
