#include "SnowTools/BreakPoint.h"

#include <getopt.h>
#include <iomanip>
#include <cassert>
#include "SnowTools/gzstream.h"

#include <boost/regex.hpp>

#define T_SPLIT_BUFF 15
#define N_SPLIT_BUFF 8
#define LOD_CUTOFF 8
#define DBCUTOFF 15
#define NODBCUTOFF 8

namespace SnowTools {

  BreakEnd::BreakEnd(const GenomicRegion& g, int mq, const std::string& chr_n) {
    gr = g;
    mapq = mq;
    chr_name = chr_n;
    cpos = -1; nm = -1; matchlen = -1;
  }
  
  ReducedBreakEnd::ReducedBreakEnd(const GenomicRegion& g, int mq, const std::string& chr_n) {
    gr = g;
    mapq = mq;
    chr_name = chr_n;
  }

  bool BreakPoint::operator==(const BreakPoint &bp) const {
    return (b1.gr == bp.b1.gr && b2.gr == bp.b2.gr && bp.insertion == insertion); 
  }
    
  std::string BreakPoint::getHashString() const {
    
    bool isdel = insertion.length() == 0;
    std::string st = std::to_string(b1.gr.chr) + "_" + std::to_string(b1.gr.pos1) + "_" + std::to_string(this->getSpan()) + (isdel ? "D" : "I");
    return st;
  }
  
  BreakEnd::BreakEnd(const BamRead& b) {
    gr.chr = b.ChrID(); 
    gr.pos1 = -1;
    gr.pos2 = -1;
    cpos = -1;
    mapq = b.MapQuality();
    chr_name = b.GetZTag("MC"); 
    assert(chr_name.length());
    assert(chr_name != "23");
    nm = std::max(b.GetIntTag("NM") - (int)b.MaxInsertionBases() - (int)b.MaxDeletionBases(), 0);
  }

  void BreakPoint::order() {

    if (b1.gr < b2.gr)
      return;
    
    std::swap(b1, b2);
    //flip(b1, b2);
    
  }

  bool BreakPoint::valid() const {

    // debug
    return true;
    
    if (!(b1.gr.strand == '+' || b1.gr.strand == '-') || !(b2.gr.strand == '+' || b2.gr.strand == '-')) {
      std::cerr << "b1.strand " << b1.gr.strand << " b2.strand " << b2.gr.strand << std::endl;
      return false;
    }

    // b1 is less than b2, or the same but signifies inverted connection
    if ((b1.gr < b2.gr) || (b1.gr.chr == b2.gr.chr && b1.gr.pos1 == b2.gr.pos1)) 
      return true;
    
    std::cerr << b1.gr << " " << b2.gr << std::endl;
    return false;
  }

  int BreakPoint::getSpan() const { 
    if (num_align == 1 && insertion == "") {// deletion
      return (abs((int)b1.gr.pos1 - (int)b2.gr.pos1) - 1);
    }
    if (num_align == 1) 
      return (insertion.length()); // insertion
    if (b1.gr.chr == b2.gr.chr)
      return abs((int)b1.gr.pos1-(int)b2.gr.pos1);
    else
      return -1;
  }

  std::string BreakEnd::hash() const {
    
    return (std::to_string(gr.chr) + ":" + std::to_string(gr.pos1));
    
  }

}
