#include "SnowTools/GenomicRegion.h"

#include <cassert>
#include "SnowTools/gzstream.h" 

// 4 billion
#define END_MAX 4000000000

namespace SnowTools {

// return the width of the genomic region
int GenomicRegion::width() const {
  return pos2 - pos1 + 1;
}

// returns 0 for no overlaps, 1 for partial and 2 for complete
int GenomicRegion::getOverlap(const GenomicRegion& gr) const {

  if (gr.chr != chr)
    return 0;
  
  // argument pos1 is in
  bool gr1_in = gr.pos1 >= pos1 && gr.pos1 <= pos2;
  // argument pos2 is in
  bool gr2_in = gr.pos2 >= pos1 && gr.pos2 <= pos2;
  // object pos1 is in
  bool pos1_in = pos1 >= gr.pos1 && pos1 <= gr.pos2;
  // object pos2 is in
  bool pos2_in = pos2 >= gr.pos1 && pos2 <= gr.pos2;

  // object is in the argument
  if (pos1_in && pos2_in) 
    return 3;

  // argument is in the oboject
  if ( gr1_in && gr2_in)
    return 2;

  // partial overlap
  if (gr1_in || gr2_in || pos1_in || pos2_in)
    return 1;

  return 0;

}


std::string GenomicRegion::ChrName(const bam_hdr_t* h) const {

  std::string cc;
  if (h) {
    if (chr >= h->n_targets)
      throw std::invalid_argument( "GenomicRegion::ChrName - not enough targets in BAM hdr to cover ref id");
    else
      cc = std::string(h->target_name[chr]);
  } else {
    cc = chrToString(chr);
  }
  return cc;
}

  std::string GenomicRegion::ChrName(const BamHeader& h) const {
    
    std::string cc;
    if (!h.isEmpty()) {
      if (chr >= h.NumSequences())
	throw std::invalid_argument( "GenomicRegion::ChrName - not enough targets in BamHeader to cover ref id");
      else
	cc = h.IDtoName(chr); // std::string(h->target_name[chr]);
    } else {
      cc = chrToString(chr);
    }
    return cc;
  }

  
// write genomic region to a string
std::string GenomicRegion::toString() const {
  std::stringstream out;
  out << chrToString(chr) << ":" << SnowTools::AddCommas<int>(pos1) << "-" << AddCommas<int>(pos2) << "(" << 
    strand << ")"; 
  return out.str();
}

  std::string GenomicRegion::pointString() const {
    std::stringstream out;
    out << chrToString(chr) << ":" << SnowTools::AddCommas<int>(pos1) << "(" << strand << ")";
    return out.str();
  }

void GenomicRegion::pad(int32_t pad) {

  if (-pad*2 > width())
    throw std::out_of_range(
         "GenomicRegion::pad - negative pad values can't obliterate GenomicRegion with val " + 
	 std::to_string(chr) + ":" + std::to_string(pos1) + "-" + std::to_string(pos2) + 
	 " and pad " + std::to_string(pad));

  pos1 -= pad;
  pos2 += pad;

  //if (pad > pos1)
  //  pos1 = 1;
  //else
  //  pos1 = pos1-pad;

  //const int32_t maxpos = 250000000;
  //pos2 = std::min(pos2+pad, maxpos); // 2500000000 is dummy for now. should be chr end

}

bool GenomicRegion::operator<(const GenomicRegion& b) const {
  return (chr < b.chr) || (chr == b.chr && pos1 < b.pos1) || (chr==b.chr && pos1 == b.pos1 && pos2 < b.pos2);
}

bool GenomicRegion::operator==(const GenomicRegion &b) const {
  return (chr == b.chr && pos1 == b.pos1 && b.pos2 == pos2);
}

bool GenomicRegion::operator<=(const GenomicRegion &b) const {
  return (*this < b || *this == b);
}

  std::string GenomicRegion::toPrettyString() const {
    
    std::stringstream ss;
    ss << (chr + 1) << ":" << AddCommas(pos1) << "-" << AddCommas(pos2);
    return ss.str();
    
  }

std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr) {
  out << gr.toString();
  return out;
}

GenomicRegion::GenomicRegion(const std::string& reg, bam_hdr_t* h) 
{
  
  if (h == nullptr)
    std::cerr <<" NULL HEADER in GenomicRegion::GenomicRegion(string, bam_hdr_t *): " << std::endl;

  // scrub String
  std::string reg2 = SnowTools::scrubString(reg, "chr");

  // use htslib region parsing code
  int tid, beg, end;
  const char * q = hts_parse_reg(reg2.c_str(), &beg, &end);
  if (q) {
    char *tmp = (char*)alloca(q - reg2.c_str() + 1); // stack alloc
    strncpy(tmp, reg2.c_str(), q - reg2.c_str());
    tmp[q - reg2.c_str()] = 0;
    tid = bam_name2id(h, tmp);
    if (tid < 0) {
      std::string inv = "GenomicRegion constructor: Failed to set region for " + reg;
      throw std::invalid_argument(inv);
    }
    
    // check that it wasn't a single region, but if so, fix (e.g. 1:1 gets fixed to 1:1-1)
    if (end == 2147483647) // at max size = not set. set to pos1
      end = beg+1;

  } else {
    std::string inv = "GenomicRegion constructor: Failed to set region for " + reg;
    throw std::invalid_argument(inv);
    //tid = bam_name2id(h, reg2.c_str());
    //beg = 0;
    //end = END_MAX;
  }
  
  chr = tid;
  pos1 = beg+1;
  pos2 = end;

}

// constructor to take a pair of coordinates to define the genomic interval
GenomicRegion::GenomicRegion(int32_t t_chr, int32_t t_pos1, int32_t t_pos2, char t_strand) {

  if (t_pos2 < t_pos1 )
    throw std::invalid_argument( "GenomicRegion constructor: end pos must be >= start pos" );

  if ( !(t_strand == '+' || t_strand == '-' || t_strand == '*') )
    throw std::invalid_argument( "GenomicRegion constructor: strand must be one of +, -, *" );

  chr = t_chr;
  pos1 = t_pos1;
  pos2 = t_pos2;
  strand = t_strand;

}

std::string GenomicRegion::chrToString(int32_t ref) {

  std::string ref_id;
  if (ref < 0)
    ref_id = std::to_string(ref);
  //throw std::invalid_argument( "GenomicRegion::chrToString - ref id must be >= 0" );

  if (ref == 22)
    ref_id = "X";
  else if (ref == 23)
    ref_id = "Y";
  else if (ref == 24)
    ref_id = "M";
  else if (ref >= 0)
    ref_id = std::to_string(ref+1);
  assert(ref_id != "23");
  return ref_id;
}

// checks whether a GenomicRegion is empty
bool GenomicRegion::isEmpty() const {
  return chr == -1 && pos1 == 0 && pos2 == 0;
}


int32_t GenomicRegion::distanceBetweenStarts(const GenomicRegion &gr) const {

  if (gr.chr != chr)
    return -1;
  else
    return std::abs(pos1 - gr.pos1);//((pos1 > gr.pos1) ? (pos1 - gr.pos1) : (gr.pos1 - pos1));

}

int32_t GenomicRegion::distanceBetweenEnds(const GenomicRegion &gr) const {

  if (gr.chr != chr)
    return -1;
  else
    return std::abs(pos2 - gr.pos2);

}


void GenomicRegion::random() {
  
  uint32_t big = rand() % SnowTools::genome_size_XY;
  //SnowTools::genRandomValue(big, SnowTools::genome_size_XY, seed);
  
  for (size_t k = 0; k < 25; ++k)
    if (big < SnowTools::CHR_CLEN[k]) {
      assert(k > 0);
      chr = --k;
      assert(big > SnowTools::CHR_CLEN[chr]);
      pos1 = big - SnowTools::CHR_CLEN[chr];
      pos2 = pos1;
      return;
    }
  
  std::cerr << "Value of " << big << " outside of expected range."  << std::endl;
  exit(EXIT_FAILURE);
  
}

  GenomicRegion::GenomicRegion(const std::string& tchr, const std::string& tpos1, const std::string& tpos2, bam_hdr_t *h)
  {
    // convert the pos strings
    try {
      pos1 = std::stoi(tpos1);
    }
    catch (...) {
      std::cerr << "GenomicRegion: error making pos1 from " << tpos1 << std::endl;
      pos1 = 0;
    }
    
    // convert the pos strings
    try {
      pos2 = std::stoi(tpos2);
    }
    catch (...) {
      std::cerr << " tchr " << tchr << " tpos1 " << tpos1 << std::endl;
      std::cerr << "GenomicRegion: error making pos2 from " << tpos2 << std::endl;
      pos2 = 0;
    }
    
      chr = -1;

      // if no header, assume that it is "standard"
      if (!h) {
	try { 
	  if (tchr == "X" || tchr == "chrX")
	    chr = 22;
	  else if (tchr == "Y" || tchr == "chrY")
	    chr = 23;
	  else 
	    chr = std::stoi(SnowTools::scrubString(tchr, "chr")) - 1;
	} catch(...) {
	  throw std::invalid_argument("GenomicRegion: error making chr from string " + tchr);
	}
	return;
      } else {
	chr = bam_name2id(h, tchr.c_str());
      }

      // TODO slow.
      //bool found = false;
      /*for (int i = 0; i < h->n_targets; ++i)
	if (strcmp(tchr.c_str(), h->target_name[i]) == 0)
	  {
	    chr = i;
	    //	    found = true;
	    break;
	  }
      */
      //debug turn this back on
      //if (!found) 
      //	std::cerr << "GenomicRegion: error, could not find matching chr in header for chr string " << tchr << std::endl;

  }
}
