#pragma once

#include <cstdint> //+11
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <string_view>

extern "C" {
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "htslib/faidx.h"

}

#include "SeqLib/SeqLibUtils.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/BamWalker.h"

constexpr char BASES[16] = {
    ' ', 'A', 'C', ' ',
    'G', ' ', ' ', ' ',
    'T', ' ', ' ', ' ',
    ' ', ' ', ' ', 'N'
};

const inline std::string cigar_delimiters = "MIDNSHPX";

constexpr uint8_t CIGTAB[255] = { 0 }; // Use constexpr if all zeroes for space and compile-time init

enum Orientation : int {
    FRORIENTATION = 0,
    FFORIENTATION = 1,
    RFORIENTATION = 2,
    RRORIENTATION = 3,
    UDORIENTATION = 4
};

namespace SeqLib {

/** Basic container for a single cigar operation
 *
 * Stores a single cigar element in a compact 32bit form (same as HTSlib).
 */
class CigarField {

  friend class Cigar;

 public:

  /** Construct the cigar op by type (MIDNSHP=X) and length 
   * @param t Cigar op (MIDNSHP=X)
   * @param l Cigar length
   * @exception Throws an invalid_argument if l <= 0 or invalid cigar op
   */
  CigarField(char opChr, uint32_t len); 

  /** Construct the cigar op from the raw sam.h uint32_t (first 4 bits op, last 28 len) */
  explicit CigarField(uint32_t f) : data(f) {}

  /** Return the raw sam.h uint32_t cigar data */
  constexpr uint32_t raw() const noexcept { return data; }

  /** Return the cigar op type (one of MIDNSHPX) as a char */
  constexpr char Type() const noexcept { return bam_cigar_opchr(data); }

  /** Return the raw sam.h uint8_t cigar type (bam_cigar_op(data)) */
  constexpr uint8_t RawType() const noexcept { return bam_cigar_op(data); }

  /** Return the length of the cigar op (eg 35M returns 35) */
  constexpr uint32_t Length() const noexcept { return bam_cigar_oplen(data); }

  /** Returns true if cigar op matches bases on the reference (MDN=X) */
  constexpr bool ConsumesReference() const noexcept {
    return (bam_cigar_type(bam_cigar_op(data)) & 2) != 0;
  }

  /** Returuns true cigar op matches bases on the query (MIS=X) */
  constexpr bool ConsumesQuery() const noexcept {
    return (bam_cigar_type(bam_cigar_op(data)) & 1) != 0;
  }  

  /** Return whether two CigarField objects have same op and len */
  constexpr bool operator==(const CigarField& other) const noexcept {
    return data == other.data;
  }  

  /** Return whether two CigarField objects have different op and/or len */
  constexpr bool operator!=(const CigarField& other) const noexcept {
    return !(*this == other);
  }

  /** Print the cigar field (eg 35M) */
  friend std::ostream& operator<<(std::ostream& out, const CigarField& c) noexcept;

 private:

  // first 4 bits hold op, last 28 hold len
  uint32_t data;
  
};

/** CIGAR for a single gapped alignment
 *
 * Constructed as a vector of CigarField objects. 
 */
 class Cigar {

 public:

   /** Construct an empty CIGAR */
   Cigar() = default;

   /** Construct from a CIGAR string 
    * @param cig CIGAR string, e.g. 54M46S
    */
   explicit Cigar(const std::string& cig);
   
   using iterator = std::vector<CigarField>::iterator;        ///< Iterator for moving between CigarField ops
   using const_iterator = std::vector<CigarField>::const_iterator; ///< Const iterator for moving between CigarField ops
   using reverse_iterator       = std::vector<CigarField>::reverse_iterator;
   using const_reverse_iterator = std::vector<CigarField>::const_reverse_iterator;
   
   iterator begin() noexcept { return m_data.begin(); }       ///< Begin iterator
   iterator end()   noexcept { return m_data.end(); }         ///< End iterator
   const_iterator begin() const noexcept { return m_data.begin(); }  ///< Const begin iterator
   const_iterator end()   const noexcept { return m_data.end(); }    ///< Const end iterator

   // Reverseiteration
   reverse_iterator rbegin()       noexcept { return m_data.rbegin(); }
   reverse_iterator rend()         noexcept { return m_data.rend(); }
   const_reverse_iterator rbegin() const noexcept { return m_data.rbegin(); }
   const_reverse_iterator rend()   const noexcept { return m_data.rend(); }
   
   /** Const reference to last cigar op */
   const CigarField& back() const noexcept { return m_data.back(); }
   /** Reference to last cigar op */
   CigarField&       back() noexcept       { return m_data.back(); }

   /** Const reference to first cigar op */
   const CigarField& front() const noexcept { return m_data.front(); }
   /** Reference to first cigar op */
   CigarField&       front() noexcept       { return m_data.front(); }

   /** Returns the number of cigar ops */
   size_t size() const noexcept { return m_data.size(); }

   /// Reserve space for at least n CigarFields to avoid reallocations
   void reserve(size_t n) { m_data.reserve(n); }
   
   /// (optional) query how many you can push without realloc
   size_t capacity() const { return m_data.capacity(); }
   
   /** Access the i'th cigar op */
   CigarField&       operator[](size_t i) noexcept       { return m_data[i]; }
   /** Access the i'th cigar op (const) */
   const CigarField& operator[](size_t i) const noexcept { return m_data[i]; }

   /** Return the number of query-consumed bases */
   int NumQueryConsumed() const noexcept;

   /** Return the number of reference-consumed bases */
   int NumReferenceConsumed() const noexcept;
   
   /** Add a new cigar op */
   void add(const CigarField& c);

   /** Return whether two Cigar objects are equivalent */
   bool operator==(const Cigar& c) const noexcept;
   
   /** Return whether two Cigar objects are not equivalent */
   bool operator!=(const Cigar& c) const { return !(c == *this); }

  /** Print cigar string (eg 35M25S) */
  friend std::ostream& operator<<(std::ostream& out, const Cigar& c) noexcept;
  
   
 private:
   
   std::vector<CigarField> m_data; // should make this simpler

 };

 typedef SeqHashMap<std::string, size_t> CigarMap;

/** Class to store and interact with a SAM alignment record
 *
 * HTSLibrary reads are stored in the bam1_t struct. Memory allocation
 * is taken care of by bam1_t init, and deallocation by destroy_bam1. This
 * class is a C++ interface that automatically takes care of memory management
 * for these C allocs/deallocs. The only member of BamRecord is a bam1_t object.
 * Alloc/dealloc is taken care of by the constructor and destructor.
 */
class BamRecord {

  friend class BWAAligner;

 public:

  /** Construct a BamRecord manually from a name, sequence, cigar and location
   * @param qname Name of the read
   * @param seq Sequence of the read (compsed of ACTG or N).
   * @param gr Location of the alignment
   * @param cig Cigar alignment
   * @exception Throws an invalid_argument exception if length of seq is not commensurate
   * with number of query-bases consumed in cigar. 
   * @exception Throws an invalid_argument exception if width of gr is not commensurate
   * with number of reference-bases consumed in cigar. 
   */
  BamRecord(std::string_view qname,
	    std::string_view seq,
	    const GenomicRegion& gr,
	    const Cigar& cig);
    
 /**
   * Take ownership of a raw bam1_t* (caller must not destroy it).
   * Wrapped in a shared_ptr so it'll be freed on destruction.
   */
  explicit BamRecord(bam1_t* raw)
    : b(raw, Bam1Deleter()) {}

  /** Check if a read is empty (not initialized)
   * @return true if read was not initialized with any values
   */
  bool isEmpty() const { return !b; }

  /** Explicitly pass a bam1_t to the BamRecord. 
   *
   * The BamRecord now controls the memory, and will delete at destruction
   * @param a An allocated bam1_t
   */
  //void assign(bam1_t* a);

  /** Make a BamRecord with no memory allocated and a null header */
  BamRecord();

  /** BamRecord is aligned on reverse strand */
  inline bool ReverseFlag() const { return b ? ((b->core.flag&BAM_FREVERSE) != 0) : false; }

  /** BamRecord has mate aligned on reverse strand */
  inline bool MateReverseFlag() const { return b ? ((b->core.flag&BAM_FMREVERSE) != 0) : false; }

  /** BamRecord has is an interchromosomal alignment */
  inline bool Interchromosomal() const { return b ? b->core.tid != b->core.mtid && PairMappedFlag() : false; }

  /** BamRecord is a duplicate */
  inline bool DuplicateFlag() const { return b ? ((b->core.flag&BAM_FDUP) != 0) : false; }

  /** BamRecord is a secondary alignment */
  inline bool SecondaryFlag() const { return b ? ((b->core.flag&BAM_FSECONDARY) != 0) : false; }

  /** BamRecord is paired */
  inline bool PairedFlag() const { return b ? ((b->core.flag&BAM_FPAIRED) != 0) : false; }

  /** Get the relative pair orientations 
   * 
   * 0 - FR (RFORIENTATION) (lower pos read is Fwd strand, higher is reverse)
   * 1 - FF (FFORIENTATION)
   * 2 - RF (RFORIENTATION)
   * 3 - RR (RRORIENTATION)
   * 4 - Undefined (UDORIENTATION) (unpaired or one/both is unmapped)
   */
  int PairOrientation() const;
  
  /** BamRecord is failed QC */
  inline bool QCFailFlag() const { return b ? ((b->core.flag&BAM_FQCFAIL) != 0) : false; }

  /** BamRecord is supplementary alignment */
  inline bool SupplementaryFlag() const { return b ? ((b->core.flag&BAM_FSUPPLEMENTARY) != 0) : false; }

  /** BamRecord is mapped */
  inline bool MappedFlag() const { return b ? ((b->core.flag&BAM_FUNMAP) == 0) : false; }

  /** BamRecord mate is mapped */
  inline bool MateMappedFlag() const { return b ? ((b->core.flag&BAM_FMUNMAP) == 0) : false; }

  /** BamRecord is mapped and mate is mapped and in pair */
  inline bool PairMappedFlag() const { return b ? (!(b->core.flag&BAM_FMUNMAP) && !(b->core.flag&BAM_FUNMAP) && (b->core.flag&BAM_FPAIRED) ) : false; }

  /** BamRecord is mapped in proper pair */
  inline bool ProperPair() const { return b ? (b->core.flag&BAM_FPROPER_PAIR) : false;} 

  /** BamRecord has proper orientation (FR) */
  bool ProperOrientation() const;

  /** Count the total number of N bases in this sequence */
  int32_t CountNBases() const;

  /** Trim the sequence down by removing bases from ends with low quality scores. Stores the
   * trimmed sequence in the GV tag, but does not affect any other part of read.
   * @param qualTrim Minimal quality score, zero-based (eg # == 2)
   * @param startpoint Returns the new starting point for the sequence
   * @param endpoint Return the new ending point for the sequence
   */
  void QualityTrimmedSequence(int32_t qualTrim, int32_t& startpoint, int32_t& endpoint) const;

  /** Retrieve the quality trimmed seqeuence from QT tag if made. Otherwise return normal seq */
  //std::string QualitySequence() const;

  /** Get the alignment position */
  inline int32_t Position() const { return b ? b->core.pos : -1; }

  /** Get the alignment position, including soft clips */
  //int32_t PositionWithSClips() const;
  
  /** Get the alignment position of mate */
  inline int32_t MatePosition() const { return b ? b->core.mpos: -1; }

  /** Count the number of secondary alignments by looking at XA tag.
   * @note A secondary alignment is an alternative mapping. This may not
   * work for non-BWA aligners that may not place the XA tag.
   */
  //int32_t CountBWASecondaryAlignments() const;

  /** Count the number of chimeric alignments by looking at XP and SA tags 
   * @note A secondary alignment is an alternative mapping. This may not
   * work for non-BWA aligners that may not place the XP/SA tags. BWA-MEM 
   * used the XP tag prior to v0.7.5, and SA aftewards.
   */
  //int32_t CountBWAChimericAlignments() const;

  /** Get the end of the alignment */
  int32_t PositionEnd() const;

  /** Get the end of the alignment, including soft clips */
  //int32_t PositionEndWithSClips() const;

  /** Get the end of the aligment mate pair */
  int32_t PositionEndMate() const;

  /** Get the chromosome ID of the read */
  int32_t ChrID() const;
  
  /** Get the chrosome ID of the mate read */
  int32_t MateChrID() const;
  
  /** Get the mapping quality */
  int32_t MapQuality() const;

  /** Set the qc fail flag on/off (true -> on) */
  void SetQCFail(bool f);
  
  /** Set the mapping quality */
  void SetMapQuality(int32_t m);

  /** Set the chr id */
  void SetChrID(int32_t i);

  /** Set the chr id of mate */
  void SetChrIDMate(int32_t i);
  
  /** Set the position of the mate read */
  void SetPositionMate(int32_t i);

  /** Set the pair mapped flag on/off (true -> on) */
  void SetPairMappedFlag(bool f);

  /** Set the mate reverse flag on/off (true -> on) */
  void SetMateReverseFlag(bool f);

  /** Get the number of cigar fields */
  inline int32_t CigarSize() const { return b ? b->core.n_cigar : -1; }
  
  /** Check if this read is first in pair */
  inline bool FirstFlag() const { return (b->core.flag&BAM_FREAD1); }

  /** Check if this read is last in pair */
  inline bool LastFlag() const { return (b->core.flag&BAM_FREAD2); }
  
  /** Get the qname of this read as a string */
  inline std::string Qname() const { return std::string(bam_get_qname(b)); }
  
  /** Get the qname of this read as a char array */
  inline char* QnameChar() const { return bam_get_qname(b); }
  
  /** Get the full alignment flag for this read */
  inline uint32_t AlignmentFlag() const { return b->core.flag; }
  
  /** Get the insert size for this read */
  inline int32_t InsertSize() const { return b->core.isize; } 

  /** Get the read group, first from qname, then by RG tag 
   * @return empty string if no readgroup found
   */
  std::string ParseReadGroup() const;

  /** Get the insert size, absolute value, and always taking into account read length */
  inline int32_t FullInsertSize() const {

    if (b->core.tid != b->core.mtid || !PairMappedFlag())
      return 0;

    return std::abs(b->core.pos - b->core.mpos) + GetCigar().NumQueryConsumed();

  }
  
  /** Get the number of query bases of this read (aka length) */
  inline int32_t Length() const { return b->core.l_qseq; }
  
  /** Append a tag with new value, delimited by 'x' */
  //void SmartAddTag(const std::string& tag, const std::string& val);
  
  /** Set the query name */
  void SetQname(std::string_view name);

  /** Set the quality scores 
   * @param quals String of quality scores or empty string
   * @param offset Offset parameter for encoding (eg 33)
   * @exception Throws an invalid_argument if n is non-empty
   * and different length than sequence
   */
  void SetQualities(std::string_view quals, int offset);

  /** Set the sequence name 
   * @param seq Sequence in upper-case (ACTGN) letters. 
   */
  void SetSequence(std::string_view seq);

  /** Set the cigar field explicitly 
   * @param c Cigar operation to set
   * @note Will not check if the cigar ops are consistent with 
   * the length of the sequence.
   */
  void SetCigar(const Cigar& c);

  /** Print a SAM-lite record for this alignment */
  friend std::ostream& operator<<(std::ostream& out, const BamRecord &r);

  /** Return read as a GenomicRegion */
  GenomicRegion AsGenomicRegion() const;

  /** Return mate read as a GenomicRegion */
  GenomicRegion AsGenomicRegionMate() const;

   /** Return the number of "aligned bases" in the same style as BamTools
    *
    * BamTools reports AlignedBases, which for example returns the literal strings (for diff CIGARs):
    * 3S5M - CTG
    * 5M - CTAGC
    * 3M1D3M - ATG-TGA 
    * 3M1I3M - ATGCTGA
    *
    * @return The number of M, D, X, = and I bases
    */
  int NumAlignedBases() const;

  /** Return the max single insertion size on this cigar */
  uint32_t MaxInsertionBases() const;

  /** Return the max single deletion size on this cigar */
  uint32_t MaxDeletionBases() const;

  /** Get the number of matched bases in this alignment */
  uint32_t NumMatchBases() const;

  /** Retrieve the CIGAR as a more managable Cigar structure */
  Cigar GetCigar() const;

  /** Retrieve the inverse of the CIGAR as a more managable Cigar structure */
  Cigar GetReverseCigar() const;

  /** Remove the sequence, quality and alignment tags. 
   * Make a more compact alignment stucture, without the string data
   */
  void ClearSeqQualAndTags();

  /** Retrieve the sequence of this read as a string (ACTGN) */
  std::string Sequence() const;

  /** Get the quality scores of this read as a string 
   * @param offset Encoding offset for phred quality scores. Default 33
   * @return Qualties scores after converting offset. If first char is empty, returns empty string
   */
  std::string Qualities(int offset = 33) const;

  /** Get the start of the alignment on the read, by removing soft-clips
   * Do this in the reverse orientation though.
   */
  int32_t AlignmentPositionReverse() const;

  /** Get the end of the alignment on the read, by removing soft-clips
   * Do this in the reverse orientation though.
   */
  int32_t AlignmentEndPositionReverse() const;

  /** Get the start of the alignment on the read, by removing soft-clips
   */
  int32_t AlignmentPosition() const;
  
  /** Get the end of the alignment on the read, by removing soft-clips
   */
  int32_t AlignmentEndPosition() const;
  
  /** Get the number of soft clipped bases */
  int32_t NumSoftClip() const;

  /** Get the number of hard clipped bases */
  int32_t NumHardClip() const;

  /** Get the number of clipped bases (hard clipped and soft clipped) */
  int32_t NumClip() const;
  
  /** Get a string (Z) tag 
   * @param tag Name of the tag. eg "XP"
   * @param out The string to be filled in with the tag information
   * @return Returns true if the tag is present, even if empty. Return false if no tag or not a Z tag.
   */
  bool GetZTag(const std::string_view tag, std::string& out) const;
  
  /** Get a string of either Z, f or i type. Useful if tag type not known at compile time.
   * @param tag Name of the tag. eg "XP"
   * @param out The string to be filled in with the tag information
   * @return Returns true if the tag is present and is either Z or i, even if empty. Return false if no tag or not Z or i.
   */  
  bool GetTag(std::string_view tag, std::string& out) const;
  
  /** Get a vector of type int from a Z tag delimited by "^"
   * Smart-tags allow one to store vectors of strings, ints or doubles in the alignment tags, and
   * do not require an additional data structure on top of bseq1_t. 
   * @param tag Name of the tag eg "AL"
   * @return A vector of ints, retrieved from the x delimited Z tag
   * @exception Throws an invalid_argument if cannot convert delimited field val to int
   */
  //std::vector<int> GetSmartIntTag(const std::string& tag) const;

  /** Get a vector of type double from a Z tag delimited by "x"
   * Smart-tags allow one to store vectors of string, ints or doubles in the alignment tags, and
   * do not require an additional data structure on top of bseq1_t. 
   * @param tag Name of the tag eg "AL"
   * @return A vector of double elems, retrieved from the "^" delimited Z tag
   * @exception Throws an invalid_argument if cannot convert delimited field val to double
   */
  //std::vector<double> GetSmartDoubleTag(const std::string& tag) const;

  /** Get a vector of strings from a Z tag delimited by "^"
   * Smart-tags allow one to store vectors of strings, ints or doubles in the alignment tags, and
   * do not require an additional data structure on top of bseq1_t. 
   * @param tag Name of the tag eg "CN"
   * @return A vector of strngs, retrieved from the x delimited Z tag
   */
  //std::vector<std::string> GetSmartStringTag(const std::string& tag) const;

  /** Get an int (i) tag 
   * @param tag Name of the tag. eg "XP"
   * @param t Value to be filled in with the tag value.
   * @return Return true if the tag exists.
   */
  bool GetIntTag(std::string_view tag, int32_t& t) const;

  /** Get a float (f) tag 
   * @param tag Name of the tag. eg "AS"
   * @param t Value to be filled in with the tag value.
   * @return Return true if the tag exists.
   */
  bool GetFloatTag(std::string_view tag, float& t) const;

  /** Add a string (Z) tag
   * @param tag Name of the tag. eg "XP"
   * @param val Value for the tag
   */
  void AddZTag(std::string_view tag, std::string_view val);

  /** Add an int (i) tag
   * @param tag Name of the tag. eg "XP"
   * @param val Value for the tag
   */
  void AddIntTag(std::string_view tag, int32_t val);

  /** Set the chr id number 
   * @param id Chromosome id. Typically is 0 for chr1, etc
   */
  void SetID(int32_t id);
  
  /** Set the alignment start position
   * @param pos Alignment start position
   */
  void SetPosition(int32_t pos);
  
  /** Convert CIGAR to a string
   */
  std::string CigarString() const;
  
  /** Return a human readable chromosome name assuming chr is indexed
   * from 0 (eg id 0 return "1")
   * @note This is a quick convienence function, and is not robust to non-numbered
   * chromosomes (eg chrX becomes 23). For accurate string representation of 
   * any chromosomes, use the full ChrName with BamHeader input.
   */
  /*  inline std::string ChrName() const {
    std::stringstream ss;
    ss << (b->core.tid + 1);

    return ss.str();
    //return std::to_string(b->core.tid + 1); //c++11
    }*/

  /** Retrieve the human readable chromosome name. 
   * @param h Dictionary for chr name lookup. If it is empty, assumes this is chr1 based reference.
   * @exception Throws an out_of_range exception if chr id is not in dictionary
   * @return Empty string if chr id < 0, otherwise chromosome name from dictionary.
   */
  std::string ChrName(const SeqLib::BamHeader& hdr) const;

  /** Return a short description (chr:pos) of this read */
  std::string Brief() const;

  /** Return a short description (chr:pos) of this read's mate */
  std::string BriefMate() const;

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

  /** Return the number of bases on the query that are covered by a match (M) on both reads 
   * This is for tracking overlapping coverage on the reads, regardless of their alignment locations.
   * For instance, two reads with 101M will have overlapping coverage of 101, regardless of alignment location.
   * A read with 50S50M and 50M50S will have 0 overlapping coverage.
   */
  int OverlappingCoverage(const BamRecord& r) const;
  
  /** Return the shared pointer */
  //SeqPointer<bam1_t> shared_pointer() const { return b; }

  // Less than operator
  bool operator<(const BamRecord& other) const;
  
  // Equality operator
  bool operator==(const BamRecord& other) const;
  
protected:
  
  SeqPointer<bam1_t> b; // bam1_t shared pointer

private:

  // init an empty bamread
  void init_();
  
};

 typedef std::vector<BamRecord> BamRecordVector; ///< Store a vector of alignment records
  
  typedef std::vector<BamRecordVector> BamRecordClusterVector; ///< Store a vector of alignment vectors
  
  /** @brief Sort methods for alignment records
   */
  namespace BamRecordSort {
    
    /** @brief Sort by read position 
     */
    struct ByReadPosition
    {
      bool operator()( const BamRecord& lx, const BamRecord& rx ) const {
	return (lx.ChrID() < rx.ChrID()) || (lx.ChrID() == rx.ChrID() && lx.Position() < rx.Position());
      }
    };
    
    /** @brief Sort by mate position 
     */
    struct ByMatePosition
    {
      bool operator()( const BamRecord& lx, const BamRecord& rx ) const {
	return (lx.MateChrID() < rx.MateChrID()) || (lx.MateChrID() == rx.MateChrID() && lx.MatePosition() < rx.MatePosition());
      }
    };
    
  }
  
}
