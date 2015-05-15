#ifndef SNOWTOOLS_GCHAIN_H__
#define SNOWTOOLS_GCHAIN_H__

#include "SnowTools/GenomicRegionCollection.h"

namespace SnowTools {

template <typename T>
class gChain {

 public:

  gChain() {}

  gChain(const GenomicRegionCollection<T>& x, const GenomicRegionCollection<T>& y);

  GenomicRegionCollection<T> lift(GenomicRegionCollection<T>& g);

  bool checkValidity() const;

 private:

  GenomicRegionCollection<T>m_galx;
  GenomicRegionCollection<T> m_galy;

  std::vector<int32_t> m_pad_left;
  std::vector<int32_t> m_pad_right;

};

}
#endif
