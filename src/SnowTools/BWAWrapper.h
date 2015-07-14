#ifndef SNOWTOOLS_BWAWRAPPER_H__
#define SNOWTOOLS_BWAWRAPPER_H__

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>

#include "SnowTools/BamRead.h"
#include "htslib/sam.h"

#define MEM_F_SOFTCLIP  0x200

extern "C" {
  #include "bwa/bwa.h"
  #include "bwa/bwt.h"
  #include "bwa/bntseq.h"
  #include "bwa/kseq.h"
  #include <stdlib.h>
  #include "bwa/utils.h"
  #include "bwa/bwamem.h"
  //#include "bwa/is.c"
  int is_bwt(ubyte_t *T, int n);

}

KSEQ_DECLARE(gzFile)

namespace SnowTools {

  /** Structure to hold unaligned sequence (name and bases)
   */
  struct USeq {
    std::string name;
    std::string seq;
  };
  typedef std::vector<USeq> USeqVector;
 
/** Calls BWA-MEM on sequence queries and returns aligned reads, all in memory 
 */
class BWAWrapper {

 public:

  /** Create a new BWAWrapper
   */
  BWAWrapper() { 
    memopt = mem_opt_init();
    memopt->flag |= MEM_F_SOFTCLIP;
  }

  /** Destroy the BWAWrapper 
   *
   * This will call the destructor on the index and options stucture
   */
  ~BWAWrapper() { 
    
    if (idx)
      bwa_idx_destroy(idx);
    if (memopt)
      free(memopt);
  }

  /** Create a bam_hdr_t from the loaded index files */
  bam_hdr_t * HeaderFromIndex() const;

  /** Convert a bns to a header string */
  std::string bwa_print_sam_hdr2(const bntseq_t *bns, const char *hdr_line) const;

  /** Construct a bam_hdr_t from a header string */
  bam_hdr_t* sam_hdr_read2(const std::string& hdr) const;

  void alignSingleSequence(const std::string& seq, const std::string& name, BamReadVector& vec, 
			   bool keep_secondary);

  /** Construct a new bwa index for this object. 
   */
  void constructIndex(const USeqVector& v);

  /** Retrieve a bwa index object from disk
   */
  void retrieveIndex(const std::string& file);

  /** Dump the stored index to files 
   */
  void writeIndexToFiles(const std::string& index_name);

  /** Return the index */
  bwaidx_t* getIndex() const { return idx; }

 private:

  mem_opt_t * memopt;

  bwaidx_t* idx = 0;

  // overwrite the bwa bwt_pac2pwt function
  bwt_t *__bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr);

  // add an anns (chr annotation structure) 
  bntann1_t* __add_to_anns(const std::string& name, const std::string& seq, bntann1_t * ann, size_t offset);

  // overwrite the bwa-mem add1 function, which takes a sequence and adds to pac
  uint8_t* __add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q);

  // make the pac structure (2-bit encoded packed sequence)
  uint8_t* __make_pac(const USeqVector& v, bool for_only, bool write_file);

  void __write_pac_to_file(const std::string& file);

  std::string __print_bns();
};

}


#endif
