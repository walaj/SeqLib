#ifndef SNOWTOOLS_READ_H__
#define SNOWTOOLS_READ_H__

#include <cstdint>

#include "SnowTools/HTSTools.h"
#include <iostream>

namespace SnowTools {

/** Class to store and interact with an HTSLib bam1_t read.
 *
 * HTSLibrary reads are stored in the bam1_t struct. Memory allocation
 * is taken care of by bam1_t init, and deallocation by destroy_bam1. This
 * class is a C++ interface that automatically takes care of memory management
 * for these C allocs/deallocs. The only member of BamRead is a bam1_t object.
 * Alloc/dealloc is taken care of by the constructor and destructor.
 */
class BamRead {

 public:
  
  /** Construct an empty BamRead by calling bam_init1() 
   */
  BamRead();

  /** Destory a BamRead by calling bam_destroy1(b)
   */
  ~BamRead();

  /** BamRead is aligned on reverse strand */
  inline bool ReverseFlag() const { return (b->core.flag&BAM_FREVERSE) != 0; }

  /** Mate of BamRead is aligned on reverse strand */
  inline bool MateReverseFlag() const { return (b->core.flag&BAM_FMREVERSE) != 0; }
  
  /** Get the alignment position */
  inline int32_t Position() const { return b->core.pos; }
  
  /** Get the alignment position of the mate */
  inline int32_t MatePosition() const { return b->core.mpos; }
  
  /** Get the chromosome ID of the read */
  inline int32_t ChrID() const { return b->core.tid; }
  
  /** Get the chrosome ID of the mate read */
  inline int32_t MateChrID() const { return b->core.mtid; }
  
  /** Check if this read is mapped*/
  inline bool MappedFlag() const { return !(b->core.flag&BAM_FUNMAP); }
  
  /** Check if this reads mate is mapped*/
  inline bool MateMappedFlag() const { return !(b->core.flag&BAM_FMUNMAP); }
  
  /** Get the mapping quality */
  inline int32_t MapQuality() const { return b->core.qual; }
  
  /** Get the number of cigar fields */
  inline int32_t CigarSize() const { return b->core.n_cigar; }
  
  /** Check if this read is first in pair */
  inline bool FirstFlag() const { return (b->core.flag&BAM_FREAD1); }
  
  /** Get the qname of this read as a string */
  inline std::string Qname() const { return std::string(bam_get_qname(b)); }
  
  /** Get the qname of this read as a char array */
  inline char* QnameChar() const { return bam_get_qname(b); }
  
  /** Get the full alignment flag for this read */
  inline int32_t AlignmentFlag() const { return b->core.flag; }
  
  /** Get the insert size for this read */
  inline int32_t InsertSize() const { return b->core.isize; } 
  
  /** Get the number of query bases of this read (aka length) */
  inline int32_t Length() const { return b->core.l_qseq; }
  
  void SetMapQuality(int32_t m) { b->core.qual = m; } 

  /** Set the query name */
  void SetQname(std::string n);

  /** Set the sequence name */
  void SetSequence(std::string s);

  friend std::ostream& operator<<(std::ostream& out, const BamRead &r);

  /** Get the sequence of this read as a string */
  inline std::string Sequence() const { 
    uint8_t * p = bam_get_seq(b);
    char c[b->core.l_qseq+1]; // need +1 for \0 at end
    for (int32_t i = 0; i < b->core.l_qseq; ++i)
      c[i] = BASES[bam_seqi(p, i)];
    c[b->core.l_qseq] = '\0'; // null terminate for proper std::string construciton
    return std::string(c);
  }
  
  /** Get the number of soft clipped bases */
  inline int32_t NumSoftClip() const {
      int32_t p = 0;
      uint32_t* c = bam_get_cigar(b);
      for (int32_t i = 0; i < b->core.n_cigar; ++i)
	if (c[i] & BAM_CSOFT_CLIP)
	  p += bam_cigar_oplen(c[i]);
      return p;
    }

  /** Get the number of clipped bases (hard clipped and soft clipped) */
  inline int32_t NumClip() const {
    int32_t p = 0;
    uint32_t* c = bam_get_cigar(b);
    for (int32_t i = 0; i < b->core.n_cigar; ++i)
      if ( (c[i] & BAM_CSOFT_CLIP) || (c[i] & BAM_CHARD_CLIP) )
	p += bam_cigar_oplen(c[i]);
    return p;
  }
  
  /** Get a string (Z) tag 
   * @param tag Name of the tag. e.g. "XP"
   * @return The value stored in the tag. Returns empty string if it does not exist.
   */
  inline std::string GetZTag(const std::string& tag) const {
    uint8_t* p = bam_aux_get(b,tag.c_str());
    if (!p)
      return "";
    return std::string(bam_aux2Z(p));
    
  }
  
  //std::string toSam(bam_hdr_t* h) const;
  
  private:

  //bam_hdr_h *hdr;
  bam1_t * b;
  

};

}
#endif
