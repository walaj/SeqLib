// BWAAligner.h
#pragma once
#include "BWAIndex.h"
#include "bwa/bwa.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/UnalignedSequence.h"

namespace SeqLib {

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
		     int maxSecondary) const;

  /// Alias for above that uses SeqLIb notation
  void alignSequence(const UnalignedSequence& us,
		     BamRecordPtrVector&           out,
		     bool                       hardclip,
		     double                     keepSecFrac,
		     int                        maxSecondary) const;
  
private:
  BWAIndexPtr   index_;
  mem_opt_t*    memopt_;
  bool          copyComment_ = false;
};
}
