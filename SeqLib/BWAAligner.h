// BWAAligner.h
#pragma once
#include "BWAIndex.h"
#include "bwa/bwa.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/UnalignedSequence.h"

namespace SeqLib {
  
  struct Alignment {
    std::string qname;
    int         rid;       // reference ID
    int64_t     pos;       // 0-based leftmost
    bool        is_rev;    // strand
    int         mapq;
    std::string cigar;
  };
  
  class BWAAligner;
  
class BWAAligner {
public:
  BWAAligner(BWAIndexPtr idx)
    : index_(std::move(idx)), memopt_(mem_opt_init())
  {
    memopt_->flag |= MEM_F_SOFTCLIP;
  }

  ~BWAAligner() {
    if (memopt_) free(memopt_);
    if(seq_buf)  free(seq_buf);
    if(core_buf) free(core_buf);
  }

  /// Set the gap-open penalty (must be >= 0)
  void SetGapOpen(int gap_open);

  /// Set the gap-extension penalty (must be >= 0)
  void SetGapExtension(int gap_ext);

  /// Set the mismatch penalty (must be >= 0)
  void SetMismatchPenalty(int mismatch);

  /// Set the Z-dropoff parameter for SW alignment (must be >= 0)
  void SetZDropoff(int zdrop);

  /// Scale all alignment scoring parameters by factor a (must be >= 0)
  void SetAScore(int a);

  /// Set the 3' clipping penalty (must be >= 0)
  void Set3primeClippingPenalty(int penalty);

  /// Set the 5' clipping penalty (must be >= 0)
  void Set5primeClippingPenalty(int penalty);

  /// Set the SW alignment bandwidth (must be >= 0)
  void SetBandwidth(int bw);

  /// Set the reseed trigger threshold (must be >= 0)
  void SetReseedTrigger(float trigger);
  
  void alignSequence(const std::string& seq,
		     const std::string& name,
		     BamRecordPtrVector& out,
		     bool hardclip,
		     double keepSecFrac,
		     int maxSecondary);

  /// Alias for above that uses SeqLIb notation
  void alignSequence(const UnalignedSequence& us,
		     BamRecordPtrVector&           out,
		     bool                       hardclip,
		     double                     keepSecFrac,
		     int                        maxSecondary);

  void allocBuffer(size_t read_len);

  void alignBatch(const UnalignedSequenceVector& inputs,
		  BamRecordPtrVector& out,
		  bool hardclip, double keepSecFrac, int maxSecondary) const;

  std::vector<Alignment> alignBatch(
						const UnalignedSequenceVector& inputs);

  
private:
  BWAIndexPtr   index_;
  mem_opt_t*    memopt_;
  bool          copyComment_ = false;

  char* seq_buf = nullptr; 
  void* core_buf = nullptr;
  size_t read_len_ = 1024;

  bam1_t* convertToBam1(const UnalignedSequence& us) const;
  
  mem_alnreg_v mem_align1_reuse(const mem_opt_t* opt,
				const bwt_t* bwt,
				const bntseq_t* bns,
				const uint8_t* pac,
				int l_seq,
				const char* seq_in,
				char* seq_buf,
				void* core_buf);

  bseq1_t convertToBseq(const UnalignedSequence& us);
  
  int64_t processedCount_ = 0;
  
};
}
