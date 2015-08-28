#ifndef SNOWTOOLS_READ_H__
#define SNOWTOOLS_READ_H__

#include <cstdint>
#include <vector>
#include <iostream>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <algorithm>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "htslib/faidx.h"

#include "SnowTools/GenomicRegion.h"

static const char BASES[16] = {' ', 'A', 'C', ' ',
                               'G', ' ', ' ', ' ', 
                               'T', ' ', ' ', ' ', 
                               ' ', ' ', ' ', 'N'};

//#include "SnowTools/HTSTools.h"

namespace SnowTools {

enum class Base { A = 1, C = 2, G = 4, T = 8, N = 15 };

/*! Basic container for cigar data. 
 */
struct CigarField {
  
  char Type; /**< The cigar operation (MIDSHPN) */
  uint32_t Length; /**< The associated length */

  CigarField(char t, uint32_t l) : Type(t), Length(l) {}

  friend std::ostream& operator<<(std::ostream& out, const CigarField& c) { out << c.Type << c.Length; return out; }
  
};


typedef std::vector<CigarField> Cigar;
typedef std::unordered_map<std::string, size_t> CigarMap;

/** Class to store and interact with an HTSLib bam1_t read.
 *
 * HTSLibrary reads are stored in the bam1_t struct. Memory allocation
 * is taken care of by bam1_t init, and deallocation by destroy_bam1. This
 * class is a C++ interface that automatically takes care of memory management
 * for these C allocs/deallocs. The only member of BamRead is a bam1_t object.
 * Alloc/dealloc is taken care of by the constructor and destructor.
 */
class BamRead {

  friend class BWAWrapper;

 public:
  
  /** Construct an empty BamRead by calling bam_init1() 
   */
  void init();

  /** Explicitly pass a bam1_t to the BamRead. 
   *
   * The BamRead now controls the memory, and will delete at destruction
   * @param a An allocated bam1_t
   */
  void assign(bam1_t* a);

  /** Make a BamRead with no memory allocated and a null header */
  BamRead() : m_hdr(nullptr) {}

  /** BamRead is aligned on reverse strand */
  inline bool ReverseFlag() const { return (b->core.flag&BAM_FREVERSE) != 0; }

  /** BamRead has mate aligned on reverse strand */
  inline bool MateReverseFlag() const { return (b->core.flag&BAM_FMREVERSE) != 0; }

  /** BamRead has is an interchromosomal alignment */
  inline bool Interchromosomal() const { return b->core.tid != b->core.mtid; }

  /** BamRead is a duplicate */
  inline bool DuplicateFlag() const { return (b->core.flag&BAM_FDUP) != 0; }

  /** BamRead is a secondary alignment */
  inline bool SecondaryFlag() const { return (b->core.flag&BAM_FSECONDARY) != 0; }

  /** BamRead is failed QC */
  inline bool QCFailFlag() const { return (b->core.flag&BAM_FQCFAIL) != 0; }

  /** BamRead is mapped */
  inline bool MappedFlag() const { return (b->core.flag&BAM_FUNMAP) == 0; }

  /** BamRead mate is mapped */
  inline bool MateMappedFlag() const { return (b->core.flag&BAM_FMUNMAP) == 0; }

  /** BamRead mate is mapped */
  inline bool PairMappedFlag() const { return !(b->core.flag&BAM_FMUNMAP) && !(b->core.flag&BAM_FUNMAP); }

  /** Count the total number of N bases in this sequence */
  int32_t CountNBases() const;

  /** Trim the sequence down by removing bases from ends with low quality scores */
  std::string QualityTrimmedSequence(int32_t qualTrim, int32_t& startpoint) const;

  /** Get the alignment position */
  inline int32_t Position() const { return b->core.pos; }
  
  /** Get the alignment position of mate */
  inline int32_t MatePosition() const { return b->core.mpos; }

  /** Count the number of secondary alignments by looking at XA and XP tags */
  int32_t CountSecondaryAlignments() const;

  /** Get the end of the alignment */
  inline int32_t PositionEnd() const { return bam_endpos(b.get()); }

  /** Get the chromosome ID of the read */
  inline int32_t ChrID() const { return b->core.tid; }
  
  /** Get the chrosome ID of the mate read */
  inline int32_t MateChrID() const { return b->core.mtid; }
  
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
  
  /** Append a tag with new value, delimited by 'x' */
  void SmartAddTag(const std::string& tag, const std::string& val);
  
  /** Set the mapping quality */
  void SetMapQuality(int32_t m) { b->core.qual = m; } 

  /** Set the query name */
  void SetQname(const std::string& n);

  /** Set the sequence name */
  void SetSequence(const std::string& seq);

  friend std::ostream& operator<<(std::ostream& out, const BamRead &r);

  /** Return read as a GenomicRegion */
  GenomicRegion asGenomicRegion() const;

  /** Get the max insertion size on this cigar */
  inline uint32_t MaxInsertionBases() const {
    uint32_t* c = bam_get_cigar(b);
    uint32_t imax = 0;
    for (size_t i = 0; i < b->core.n_cigar; i++) 
      if (bam_cigar_opchr(c[i]) == 'I')
	imax = std::max(bam_cigar_oplen(c[i]), imax);
    return imax;
  }

  /** Get the max deletion size on this cigar */
  inline uint32_t MaxDeletionBases() const {
    uint32_t* c = bam_get_cigar(b);
    uint32_t dmax = 0;
    for (size_t i = 0; i < b->core.n_cigar; i++) 
      if (bam_cigar_opchr(c[i]) == 'D')
	dmax = std::max(bam_cigar_oplen(c[i]), dmax);
    return dmax;
  }

  /** Get the max deletion size on this cigar */
  inline uint32_t NumMatchBases() const {
    uint32_t* c = bam_get_cigar(b);
    uint32_t dmax = 0;
    for (size_t i = 0; i < b->core.n_cigar; i++) 
      if (bam_cigar_opchr(c[i]) == 'M')
	dmax += bam_cigar_oplen(c[i]);
    return dmax;
  }


  /** Retrieve the CIGAR as a more managable Cigar structure */
  Cigar GetCigar() const {
    uint32_t* c = bam_get_cigar(b);
    Cigar cig;
    for (int k = 0; k < b->core.n_cigar; ++k) 
      cig.push_back(CigarField("MIDSSHP=XB"[c[k]&BAM_CIGAR_MASK], bam_cigar_oplen(c[k])));
    return cig;
  }

  /** Retrieve the inverse of the CIGAR as a more managable Cigar structure */
  Cigar GetReverseCigar() const {
    uint32_t* c = bam_get_cigar(b);
    Cigar cig;
    for (int k = b->core.n_cigar - 1; k >= 0; --k) 
      cig.push_back(CigarField("MIDSSHP=XB"[c[k]&BAM_CIGAR_MASK], bam_cigar_oplen(c[k])));
    return cig;
  }

  /** Remove the sequence, quality and alignment tags */
  void clearSeqQualAndTags();

  /** Get the sequence of this read as a string */
  inline std::string Sequence() const { 
    uint8_t * p = bam_get_seq(b);
    std::stringstream ss;
    for (int32_t i = 0; i < b->core.l_qseq; ++i) 
      ss << BASES[bam_seqi(p, i)];
    return ss.str();
  }

  /** Get the quality scores of this read as a string */
  inline std::string Qualities() const { 
    uint8_t * p = bam_get_qual(b);
    std::stringstream ss;
    for (int32_t i = 0; i < b->core.l_qseq; ++i) 
      ss << (char)(p[i]+ 33);
    return ss.str();
  }

  /** Get the start of the alignment on the read, by removing soft-clips
   * Do this in the reverse orientation though.
   */
  inline int32_t AlignmentPositionReverse() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t p = 0;
    for (int32_t i = b->core.n_cigar - 1; i >= 0; --i) {
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H'))
	p += bam_cigar_oplen(c[i]);
      else // not a clip, so stop counting
	break;
    }
    return p;
  }
  
  /** Get the end of the alignment on the read, by removing soft-clips
   * Do this in the reverse orientation though.
   */
  inline int32_t AlignmentEndPositionReverse() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t p = 0;
    for (int32_t i = 0; i < b->core.n_cigar; ++i) { // loop from the end
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H'))
	p += bam_cigar_oplen(c[i]);
      else // not a clip, so stop counting
	break;
    }
    return (b->core.l_qseq - p);
  }


  /** Get the start of the alignment on the read, by removing soft-clips
   */
  inline int32_t AlignmentPosition() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t p = 0;
    for (int32_t i = 0; i < b->core.n_cigar; ++i) {
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H'))
	p += bam_cigar_oplen(c[i]);
      else // not a clip, so stop counting
	break;
    }
    return p;
  }
  
  /** Get the end of the alignment on the read, by removing soft-clips
   */
  inline int32_t AlignmentEndPosition() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t p = 0;
    for (int32_t i = b->core.n_cigar - 1; i >= 0; --i) { // loop from the end
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H'))
	p += bam_cigar_oplen(c[i]);
      else // not a clip, so stop counting
	break;
    }
    return (b->core.l_qseq - p);
  }

  /** Get the number of soft clipped bases */
  inline int32_t NumSoftClip() const {
      int32_t p = 0;
      uint32_t* c = bam_get_cigar(b);
      for (int32_t i = 0; i < b->core.n_cigar; ++i)
	if (bam_cigar_opchr(c[i]) == 'S')
	  p += bam_cigar_oplen(c[i]);
      return p;
    }

  /** Get the number of hard clipped bases */
  inline int32_t NumHardClip() const {
      int32_t p = 0;
      uint32_t* c = bam_get_cigar(b);
      for (int32_t i = 0; i < b->core.n_cigar; ++i) 
	if (bam_cigar_opchr(c[i]) == 'H')
	  p += bam_cigar_oplen(c[i]);
      return p;
    }


  /** Get the number of clipped bases (hard clipped and soft clipped) */
  inline int32_t NumClip() const {
    int32_t p = 0;
    uint32_t* c = bam_get_cigar(b);
    for (int32_t i = 0; i < b->core.n_cigar; ++i)
      if ( (bam_cigar_opchr(c[i]) == 'S') || (bam_cigar_opchr(c[i]) == 'H') )
	p += bam_cigar_oplen(c[i]);
    return p;
  }
  
  /** Get a string (Z) tag 
   * @param tag Name of the tag. e.g. "XP"
   * @return The value stored in the tag. Returns empty string if it does not exist.
   */
  inline std::string GetZTag(const std::string& tag) const {
    uint8_t* p = bam_aux_get(b.get(),tag.c_str());
    if (!p)
      return "";
    return std::string(bam_aux2Z(p));
  }
  
  /** Get a vector of ints from a Z tag delimited by "x"
   * @param tag Name of the tag e.g. "AL"
   * @return A vector of ints, retrieved from the x delimited Z tag
   */
  std::vector<int> GetSmartIntTag(const std::string& tag) const;

  /** Get a vector of strings from a Z tag delimited by "x"
   * @param tag Name of the tag e.g. "CN"
   * @return A vector of strngs, retrieved from the x delimited Z tag
   */
  std::vector<std::string> GetSmartStringTag(const std::string& tag) const;

  /** Get an int (i) tag 
   * @param tag Name of the tag. e.g. "XP"
   * @return The value stored in the tag. Returns 0 if it does not exist.
   */
  inline int32_t GetIntTag(const std::string& tag) const {
    uint8_t* p = bam_aux_get(b.get(),tag.c_str());
    if (!p)
      return 0;
    return bam_aux2i(p);
  }


  /** Add a string (Z) tag
   * @param tag Name of the tag. e.g. "XP"
   * @param val Value for the tag
   */
  inline void AddZTag(std::string tag, std::string val) {
    bam_aux_append(b.get(), tag.data(), 'Z', val.length()+1, (uint8_t*)val.c_str());
  }

  /** Add an int (i) tag
   * @param tag Name of the tag. e.g. "XP"
   * @param val Value for the tag
   */
  inline void AddIntTag(const std::string& tag, int32_t val) {
    bam_aux_append(b.get(), tag.data(), 'i', 4, (uint8_t*)&val);
  }

  /** Set the chr id number 
   * @param id Chromosome id. Typically is 0 for chr1, etc
   */
  inline void SetID(int32_t id) {
    b->core.tid = id;
  }
  
  /** Set the alignment start position
   * @param pos Alignment start position
   */
  inline void SetPosition(int32_t pos) {
    b->core.pos = pos;
  }

  /** Convert CIGAR to a string
   */
  inline std::string CigarString() const {
    std::stringstream cig;
    uint32_t* c = bam_get_cigar(b);
    for (int k = 0; k < b->core.n_cigar; ++k)
      cig << bam_cigar_oplen(c[k]) << "MIDNSHP=XB"[c[k]&BAM_CIGAR_MASK];
    return cig.str();
  }
  
  /** Retrieve the human readable chromosome name. 
   * 
   * Note that this requires that the header not be empty. If
   * it is empty, assumes this ia chr1 based reference
   */
  inline std::string ChrName() const {

    // if we have the header, convert
    if (m_hdr) {
      if (b->core.tid < m_hdr->n_targets)
	return std::string(m_hdr->target_name[b->core.tid]);
      else
	return "CHR_ERROR";
    }

    // no header, assume zero based
    return std::to_string(b->core.tid + 1);
    
  }

  /** Return a short description (chr:pos) of this read */
  inline std::string Brief(bam_hdr_t * h = nullptr) const {
    if (!h)
      return(std::to_string(b->core.tid + 1) + ":" + AddCommas<int32_t>(b->core.pos) + "(" + ((b->core.flag&BAM_FREVERSE) != 0 ? "+" : "-") + ")");
    else
      return(std::string(h->target_name[b->core.tid]) + ":" + AddCommas<int32_t>(b->core.pos) + "(" + ((b->core.flag&BAM_FREVERSE) != 0 ? "+" : "-") + ")");      
  }

  /** Return a short description (chr:pos) of this read's mate */
  inline std::string BriefMate(bam_hdr_t * h = nullptr) const {
    if (!h)
      return(std::to_string(b->core.mtid + 1) + ":" + AddCommas<int32_t>(b->core.mpos) + "(" + ((b->core.flag&BAM_FMREVERSE) != 0 ? "+" : "-") + ")");
    else
      return(std::string(h->target_name[b->core.mtid]) + ":" + AddCommas<int32_t>(b->core.mpos) + "(" + ((b->core.flag&BAM_FMREVERSE) != 0 ? "+" : "-") + ")");      
  }


  /** Strip a particular alignment tag 
   * @param tag Tag to remove
   */
  inline void RemoveTag(const char* tag) {
    uint8_t* p = bam_aux_get(b.get(), tag);
    if (p)
      bam_aux_del(b.get(), p);
  }

  /** Strip all of the alignment tags */
  inline void RemoveAllTags() {
    size_t keep = (b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;
    b->data = (uint8_t*)realloc(b->data, keep); // free the end, which has aux data
    b->l_data = keep;
    b->m_data = b->l_data;
  }

  /** Return the raw pointer */
  inline bam1_t* raw() const { return b.get(); }
  
  //std::string toSam(bam_hdr_t* h) const;
  
  private:

  //bam_hdr_h *hdr;
  //bam1_t * b;
  std::shared_ptr<bam1_t> b;
  std::shared_ptr<bam_hdr_t> m_hdr;
};

 typedef std::vector<BamRead> BamReadVector; 
 
 typedef std::vector<BamReadVector> BamReadClusterVector;

 namespace BamReadSort {

   // sort by read position
   struct ByReadPosition
   {
     bool operator()( const BamRead& lx, const BamRead& rx ) const {
       return (lx.ChrID() < rx.ChrID()) || (lx.ChrID() == rx.ChrID() && lx.Position() < rx.Position());
     }
   };

   // sort by read mate position
   struct ByMatePosition
   {
     bool operator()( const BamRead& lx, const BamRead& rx ) const {
       return (lx.MateChrID() < rx.MateChrID()) || (lx.MateChrID() == rx.MateChrID() && lx.MatePosition() < rx.MatePosition());
     }
   };

}

}
#endif
