#include "SnowTools/gChain.h"

#include <cassert>

using namespace SnowTools;

template <typename T>
gChain<T>::gChain(const GenomicRegionCollection<T>& x, const GenomicRegionCollection<T>& y) : m_galx(x), m_galy(y)
{

  m_pad_left = std::vector<int32_t>(x.size(), 0);
  m_pad_right = std::vector<int32_t>(x.size(), 0);

  if (!checkValidity())
    throw 20;

}

template <typename T>
GenomicRegionCollection<T> gChain<T>::lift(GenomicRegionCollection<T>& g)
{

  //RCODE: hits <- gr.findoverlaps(x, .Object@.galx, verbose = verbose, pintersect=pval, ...) # pairs of matches
  std::vector<size_t> query_id, subject_id;
  g.createTreeMap();
  GenomicRegionCollection<T> hits = g.findOverlaps(m_galx, query_id, subject_id);
  assert(query_id.size() == subject_id.size());

  // debug
  for (auto& i : subject_id)
    std::cout << " subjet id " << i << std::endl;

  //debug
  for (auto& i : hits)
    std::cout << " lift hits " << i << std::endl;

  //RCODE: pad.left.x = .Object@.pad.left[values(hits)$subject.id]
  std::vector<int32_t> pad_left_y(hits.size(), 0);
  std::vector<int32_t> pad_right_y(hits.size(), 0);

  //RCODE: qr.hits = .Object@.galy[values(hits)$subject.id];
  GenomicRegionCollection<T> qr_hits;
  for (auto& i : subject_id) 
    qr_hits.add(m_galy.at(i));
  
  //debug
  for (auto& i : m_galx)
    std::cout << "GALX " << i << std::endl;

  //RCODE: link.starts = start(.Object@.galx)[values(hits)$subject.id]
  //RCODE: link.ends = end(.Object@.galx)[values(hits)$subject.id]
  std::vector<int32_t> link_starts(hits.size(), 0);
  std::vector<int32_t> link_ends(hits.size(), 0);
  for (size_t i = 0; i < hits.size(); ++i) 
    {
      link_starts[i] = m_galx.at(subject_id[i]).pos1;
      link_ends[i]   = m_galx.at(subject_id[i]).pos2;
    }

  //debug
  for (size_t i = 0; i < hits.size(); ++i)
    std::cout << " link_starts " << link_starts[i] << " ends " << link_ends[i] << std::endl;

  //RCODE:  qd.overlap = ranges(hits)
  //   s.overlap is vec of 1 for scale=1, length = length(hits)
  //RCODE: shift1 = (start(qd.overlap) - link.starts)*abs(s.overlap)
  //RCODE: shift2 = ((end(qd.overlap) - link.starts + 1)*(abs(s.overlap)))-1
  std::vector<int32_t> shift1, shift2;
  assert(hits.size() == link_starts.size());
  for (size_t i = 0; i < hits.size(); ++i) 
    {
      shift1.push_back((int32_t)hits.at(i).pos1 - link_starts[i]);
      shift2.push_back((int32_t)hits.at(i).pos2 - link_starts[i]);
    }

  //debug
  for (size_t i = 0; i < hits.size(); ++i)
    std::cout << " shift1 " << shift1[i] << " shift2 " << shift2[i] << std::endl;

  std::vector<int32_t> starts(hits.size(), 0); 
  std::vector<int32_t> ends(hits.size(), 0);
  for (size_t i = 0; i < hits.size(); ++i)
    {
      if (hits.at(i).strand)
	{
	  //RCODE: starts[!neg.map] <- start(qr.hits)[!neg.map] + (shift1[!neg.map] - pad.left.y[!neg.map]);
	  starts[i] = qr_hits.at(i).pos1 + shift1[i] - pad_left_y[i];
	  //RCODE: ends[!neg.map]   <- start(qr.hits)[!neg.map] + (shift2[!neg.map] - pad.left.y[!neg.map]);
	  ends[i] = qr_hits.at(i).pos1 + shift2[i] - pad_left_y[i];	
	}
      else
	{
	  //RCODE: starts[neg.map] = end(qr.hits)[neg.map] - (shift2[neg.map] - pad.right.y[neg.map])
	  starts[i] = qr_hits.at(i).pos2 - (shift2[i] - pad_right_y[i]);
	  //RCODE: ends[neg.map] = end(qr.hits)[neg.map] - (shift1[neg.map] - pad.right.y[neg.map])
	  ends[i] = qr_hits.at(i).pos2 - (shift1[i] - pad_right_y[i]);
	}
    }

  //debug
  for (size_t i = 0; i < hits.size(); ++i)
    std::cout << " starts " << starts[i] << " ends " << ends[i] << std::endl;
      
  //RCODE: out = GRanges(seqnames(qr.hits), IRanges(starts, ends), strand = strand(qr.hits), seqlengths = seqlengths(qr.hits));  
  GenomicRegionCollection<T> out;
  for (size_t i = 0; i < hits.size(); ++i)
    {
      GenomicRegion gr(qr_hits.at(i).chr, starts[i], ends[i], qr_hits.at(i).strand);
      out.add(gr);
    }
      
  return out;
}

template <typename T>
bool gChain<T>::checkValidity() const
{
  
  if (m_galx.size() != m_galy.size())
    {
      std::cerr << "The two GenomicRegionCollection elements must have the same number of regions" << std::endl;
      return false;
    }

  for (size_t i = 0; i < m_galx.size(); ++i)
    {
      if (m_galx.at(i).width() != m_galy.at(i).width())
	{
	  std::cerr << "Each pair-wise element of the regions must have the same width" << std::endl;
	  return false;
	}
    }
  
  return true;

}
