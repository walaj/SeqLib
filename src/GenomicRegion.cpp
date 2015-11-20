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
int GenomicRegion::getOverlap(const GenomicRegion gr) const {

  if (gr.chr != chr)
    return 0;
  
  bool gr1_in = gr.pos1 >= pos1 && gr.pos1 <= pos2;
  bool gr2_in = gr.pos2 >= pos1 && gr.pos2 <= pos2;
  bool pos1_in = pos1 >= gr.pos1 && pos1 <= gr.pos2;
  bool pos2_in = pos2 >= gr.pos1 && pos2 <= gr.pos2;

  if ( (gr1_in && gr2_in) || (pos1_in && pos2_in) )
    return 2;

  if (gr1_in || gr2_in || pos1_in || pos2_in)
    return 1;

  return 0;

}


  std::string GenomicRegion::ChrName(const bam_hdr_t* h) const {
    std::string cc;
    if (h) {
      if (chr >= h->n_targets)
	std::cerr << "chr " << chr << " is bigger than provided targets of " << h->n_targets << std::endl;
      else
	cc = std::string(h->target_name[chr]);
    } else {
      cc = chrToString(chr);
    }
    return cc;
  }
  
// write genomic region to a string
std::string GenomicRegion::toString() const {
  std::stringstream out;
  //out << chrToString(chr)  << ":" << SnowUtils::AddCommas<int>(pos1) << "-" << SnowUtils::AddCommas<int>(pos2) << "(" << strand << ")"; 
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
  if (pad > pos1)
    pos1 = 1;
  else
    pos1 = pos1-pad;

  const int32_t maxpos = 250000000;
  pos2 = std::min(pos2+pad, maxpos); // 2500000000 is dummy for now. should be chr end
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

std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr) {
  out << gr.toString();
  return out;
}

GenomicRegion::GenomicRegion(const std::string& reg, bam_hdr_t* h) 
{
  // scrubtString
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
      std::cerr << "Failed to set region for region string " << reg << std::endl;
      std::cerr << "chr-id " << tid << " pos1 " << beg << " pos2 " << end << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    tid = bam_name2id(h, reg2.c_str());
    beg = 0;
    end = END_MAX;
    std::cerr << "LTID " << tid << std::endl;
  }

  chr = tid;
  pos1 = beg+1;
  pos2 = end;

}

// constructor for SnowTools::GenomicRegion that takes strings. Assumes chr string is in 
// natural (1, ..., X) or (chr1, ..., chrX) format. That is, it converts to
// BamTools format with a -1 operation.
/*GenomicRegion::GenomicRegion(std::string t_chr, std::string t_pos1, std::string t_pos2) {

  chr = GenomicRegion::chrToNumber(t_chr);
  try {
    t_pos1 = SnowTools::scrubString(t_pos1, ",");
    t_pos2 = SnowTools::scrubString(t_pos2, ",");
    pos1 = std::stoi(t_pos1);
    pos2 = std::stoi(t_pos2);
  } catch (...) { 
    std::cerr << "stoi failed in GenomicRegion constructor. Tried: " << t_pos1 << " " << t_pos2 << std::endl;
  }
}
*/

// constructor to take a pair of coordinates to define the genomic interval
GenomicRegion::GenomicRegion(int32_t t_chr, uint32_t t_pos1, uint32_t t_pos2, char t_strand) {
  chr = t_chr;
  pos1 = t_pos1;
  pos2 = t_pos2;
  strand = t_strand;
}

// convert a chromosome string into a number
int GenomicRegion::chrToNumber(std::string ref) {

  // remove the chr identifier if it is there
  if (ref.find("chr") != std::string::npos)
    ref = ref.substr(3, ref.size() - 3);

  std::string ref_id = ref;
  if (ref_id == "X")
    ref_id = "23";
  else if (ref_id == "Y")
    ref_id = "24";
  else if (ref_id == "M" || ref_id == "MT")
    ref_id = "25";
  
  int out = -1;
  try {
    out = std::stoi(ref_id);
  } catch (...) {
    //cerr << "Caught error trying to convert " << ref << " to number" << endl;
  }

  //assert(out > 0);
  return (out-1); // offset by one becuase chr1 = 0 in BamAlignment coords
}

// convert a chromosome number to a string. Assumes 
// a natural ordering (1, ...), not BamTools ordering (0, ...)
std::string GenomicRegion::chrToString(int32_t ref) {
  std::string ref_id;
  if (ref == 22)
    ref_id = "X";
  else if (ref == 23)
    ref_id = "Y";
  else if (ref == 24)
    ref_id = "M";
  else
    ref_id = std::to_string(ref+1);
  assert(ref_id != "23");
  return ref_id;
}

// checks whether a GenomicRegion is empty
bool GenomicRegion::isEmpty() const {
  return chr == -1 && pos1 == 0 && pos2 == 0;
}


//
/*uint32_t GenomicRegion::posToBigPos(int refid, int pos) {
  
  if (refid < 25)
    return 0;
  
  return CHR_CLEN[refid] + pos;
  
  }*/


int GenomicRegion::distance(const GenomicRegion &gr) const {

  if (gr.chr != chr)
    return -1;
  else
    return std::abs(pos1 - gr.pos1);//((pos1 > gr.pos1) ? (pos1 - gr.pos1) : (gr.pos1 - pos1));

}

void GenomicRegion::random(uint32_t seed) {
  
  uint32_t big;
  SnowTools::genRandomValue(big, SnowTools::genome_size_XY, seed);
  
  for (size_t k = 0; k < 25; k++)
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
	  std::cerr << "GenomicRegion: error making chr from string " << tchr << std::endl;
	}
	return;
      }

      // TODO slow.
      //bool found = false;
      for (int i = 0; i < h->n_targets; ++i)
	if (strcmp(tchr.c_str(), h->target_name[i]) == 0)
	  {
	    chr = i;
	    //	    found = true;
	    break;
	  }

      //debug turn this back on
      //if (!found) 
      //	std::cerr << "GenomicRegion: error, could not find matching chr in header for chr string " << tchr << std::endl;

  }
}
