#ifndef SNOWTOOLS_ALIGNED_CONTIG_H__
#define SNOWTOOLS_ALIGNED_CONTIG_H__

#include <algorithm>

#include "SnowTools/BamRead.h"
#include "SnowTools/BreakPoint.h"
#include "SnowTools/BamWalker.h"

namespace SnowTools {

  class AlignedContig;
  
  /*! This class contains a single alignment fragment from a contig to
   * the reference. For a multi-part mapping of a contig to the reference,
   * an object of this class represents just a single fragment from that alignment.
   */
  struct AlignmentFragment {
    
    friend AlignedContig;
    
    /*! Construct an AlignmentFragment from a BWA alignment
     * @param const reference to an aligned sequencing read
     * @param flip If the contig sequence was flipped (rev of BAM record), need to track this. This flipping occurs in AlignedContig::AlignedContig
     */
    AlignmentFragment(const BamRead &talign, bool flip, const std::unordered_set<std::string>& prefixes);
    
    //! sort AlignmentFragment objects by start position
    bool operator < (const AlignmentFragment& str) const { return (start < str.start); }

    //! print the AlignmentFragment
    std::string print() const;

    BreakEnd makeBreakEnd(bool left);
    
    /*! @function
     * @abstract Parse an alignment frag for a breakpoint
     * @param reference to a BreakPoint to be created
     * @return boolean informing whether there was a remaining indel break
     */
    bool parseIndelBreak(BreakPoint &bp);
    
    const std::vector<BreakPoint>& getIndelBreaks() const { return m_indel_breaks; }
    
    /*! Write the alignment record to a BAM file
     */
    void writeToBAM(BamWalker& bw); 

    private:

    std::vector<AlignmentFragment> secondaries;

    BamRead m_align; /**< BWA alignment to reference */

    int sub_n = 0; // number of sub optimal alignments
    
    std::vector<BreakPoint> m_indel_breaks; /**< indel variants on this alignment */
    
    Cigar m_cigar; /**< cigar oriented to assembled orientation */
    
    size_t idx = 0; // index of the cigar where the last indel was taken from 
    
    int break1 = -1; // 0-based breakpoint 1 on contig 
    int break2 = -1; /**< 0-based breakpoint 2 on contig */
    int gbreak1 = -1; /**< 0-based breakpoint 1 on reference chr */
    int gbreak2 = -1; /**< 0-based breakpoint 1 on reference chr */
    
    int start; /**< the start position of this alignment on the reference. */
    
    bool local = false; /**< boolean to note whether this fragment aligns to same location is was assembled from */
    
    AlignedContig * c; // link to the parent aligned contigs

    int di_count = 0; // number of indels

    int num_align = 0;
  };
  
  //! vector of AlignmentFragment objects
  typedef std::vector<AlignmentFragment> AlignmentFragmentVector;
  
  /*! Contains the mapping of an aligned contig to the reference genome,
   * along with pointer to all of the reads aligned to this contig, and a 
   * store of all of the breakpoints associated with this contig
   */
  class AlignedContig {
    
    friend class AlignmentFragment;
    
  public:  
    
    AlignedContig() {}
    
    AlignedContig(const BamReadVector& bav);
    
    /*! @function Determine if this contig has identical breaks and is better than another.
     * @param const reference to another AlignedContig
     * @return bool bool returning true iff this contig has identical info has better MAPQ, or equal MAPQ but longer */
    bool isWorse(const AlignedContig &ac) const;

    //! Loop through fragments and check if they overlap with window (and set local flag). Return TRUE if local found
    bool checkLocal(const GenomicRegion& window);
    
    //! return the name of the contig
    std::string getContigName() const { 
      if (!m_frag_v.size()) 
	return "";  
      return m_frag_v[0].m_align.Qname(); 
    }
    
    /*! @function get the maximum mapping quality from all alignments
     * @return int max mapq
     */
    int getMaxMapq() const { 
      int m = -1;
      for (auto& i : m_frag_v)
	if (i.m_align.MapQuality() > m)
	  m = i.m_align.MapQuality();
      return m;
      
    }
    
    /*! @function get the minimum mapping quality from all alignments
     * @return int min mapq
     */
    int getMinMapq() const { 
      int m = 1000;
      for (auto& i : m_frag_v)
	if (i.m_align.MapQuality() < m)
	  m = i.m_align.MapQuality();
      return m;
    }
    
    /*! @function dump the contigs to a fasta
     * @param ostream to write to
     */
    void printContigFasta(std::ofstream &os) const;

  /*! @function set the breakpoints on the reference by combining multi-mapped contigs
   */
  void setMultiMapBreakPairs();

  /**
   */
  void alignReads(BamReadVector &bav);

  //! return the contig sequence as it came off the assembler
  std::string getSequence() const { assert(m_seq.length()); return m_seq; }

  //! detemine if the contig contains a subsequence
  bool hasSubSequence(const std::string& subseq) const { 
    return (m_seq.find(subseq) != std::string::npos);
  }
  
  //! print this contig
  std::string print() const;

  /*! @function query if this contig contains a potential variant (indel or multi-map)
   * @return true if there is multimapping or an indel
   */
  bool hasVariant() const;

  /*! Write all of the alignment records to a BAM file
   * @param bw BamWalker opened with OpenWriteBam
   */
  void writeToBAM(BamWalker& bw);

  /*! Write all of the sequencing reads as aligned to contig to a BAM file
   * @param bw BamWalker opened with OpenWriteBam
   */
  void writeAlignedReadsToBAM(BamWalker& bw); 

  /*! @function retrieves all of the breakpoints by combining indels with global mutli-map break
   * @return vector of ind
   */
  std::vector<BreakPoint> getAllBreakPoints(bool local_restrict = true) const;

  std::vector<BreakPoint> getAllBreakPointsSecondary() const;

  std::vector<const BreakPoint*> getAllBreakPointPointers() const ;

  std::pair<int, int> getCoverageAtPosition(int pos) const;

  BamReadVector m_bamreads; // store all of the reads aligned to contig

  std::unordered_map<std::string, std::vector<int>> cov;

  std::vector<int> tum_cov, norm_cov;

  std::unordered_set<std::string> prefixes; // store the sample ids. Needed to create accurate BreakPoint genotypes

  AlignmentFragmentVector m_frag_v; // store all of the individual alignment fragments 

 private:

  std::vector<BreakPoint> m_local_breaks; // store all of the multi-map BreakPoints for this contigs 

  std::vector<BreakPoint> m_local_breaks_secondaries; // store all of the multi-map BreakPoints for this contigs 

  BreakPoint m_global_bp;  // store the single spanning BreakPoing for this contig

  std::vector<BreakPoint> m_global_bp_secondaries;  // store the single spanning BreakPoing for this contig e

  std::string m_seq = ""; // sequence of contig as it came off of assembler

};

struct PlottedRead {

  int pos;
  std::string seq;
  std::string info;

  bool operator<(const PlottedRead& pr) const {
    return (pos < pr.pos);
  }

};

typedef std::vector<PlottedRead> PlottedReadVector;

struct PlottedReadLine {

  std::vector<PlottedRead*> read_vec;
  int available = 0;
  int contig_len = 0;

  void addRead(PlottedRead *r) {
    read_vec.push_back(r);
    available = r->pos + r->seq.length() + 5;
  }

  bool readFits(PlottedRead &r) {
    return (r.pos >= available);
  }

  friend std::ostream& operator<<(std::ostream& out, const PlottedReadLine &r) {
    int last_loc = 0;
    for (auto& i : r.read_vec) {
      assert(i->pos - last_loc >= 0);
      out << std::string(i->pos - last_loc, ' ') << i->seq;
      last_loc = i->pos + i->seq.length();
    }
    int name_buff = r.contig_len - last_loc;
    assert(name_buff < 1e6);
    out << std::string(std::max(name_buff, 5), ' ');
    for (auto& i : r.read_vec) { // add the data
      out << i->info << ",";
    }
    return out;
  }

};

typedef std::vector<PlottedReadLine> PlottedReadLineVector;



}

#endif
