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

void BamReadGroup::addRead(BamRead &r)
{
  int mapqr = r.MapQuality();
  if (mapqr >=0 && mapqr <= 100)
    mapq.addElem(mapqr);
  
  int32_t this_nm = r.GetIntTag("NM");;
  //r_get_int32_tag(r, "NM", this_nm);
  if (this_nm <= 100)
    nm.addElem(this_nm);
  
  int32_t isizer = -1;
  if (!r.PairMappedFlag())
    isizer = -2;
  else if (!r.Interchromosomal())
    isizer = std::abs(r.InsertSize());
  isize.addElem(isizer);

  int32_t c = r.NumClip();
  //r_get_clip(r,c);
  clip.addElem(c);
  
  len.addElem(r.Length());

}

void BamStats::addRead(BamRead &r)
{

  // get the read group
  std::string rg = r.GetZTag("RG");
  //r_get_Z_tag(r, "RG", rg);

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
