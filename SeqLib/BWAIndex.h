// BWAIndex.h
#pragma once

#include <string>
#include <ostream>
#include <memory>

#include "SeqLib/UnalignedSequence.h"
#include "SeqLib/BamHeader.h"

extern "C" {
  #include "bwa/bwa.h"
  #include "bwa/bwt.h"
  #include "bwa/bntseq.h"
  #include "bwa/kseq.h"
  #include <stdlib.h>
  #include "bwa/utils.h"
  #include "bwa/bwamem.h"
  int is_bwt(ubyte_t *T, int n);
  KSEQ_DECLARE(gzFile)
}

namespace SeqLib {
  
/// Manages a BWA on-disk index (.bwt/.sa/.bns/.pac).
/// You can load an existing prefix or write out the in-memory index.
class BWAIndex {
public:
  BWAIndex() = default;
  ~BWAIndex();

  bool IsEmpty() const noexcept;
  
  BamHeader HeaderFromIndex() const;

  int NumSequences() const;

  std::string ChrIDToName(int id) const;

  std::string printSamHeader() const;

  void ConstructIndex(const UnalignedSequenceVector& refs);
  
  /// Load a BWA index from disk (prefix without extension).
  /// Returns true on success.
  void LoadIndex(const std::string& prefix);

  /// Write the currently loaded index back out under `prefix`.
  /// Throws std::runtime_error if no index is loaded.
  void WriteIndex(const std::string& prefix) const;

  /// Pretty-print some summary info about the loaded index.
  friend std::ostream& operator<<(std::ostream& os, const BWAIndex& idx);

  /// Pretty-print some summary info about the loaded index.
  friend std::ostream& operator<<(std::ostream& os, const BWAIndex& idx);
  
private:
  bwaidx_t* idx_ = nullptr;  ///< opaque BWA index

  /// Dump the `.pac` file for idx_->pac under `<file>.pac`
  void seqlib_write_pac_to_file(const std::string& file) const;
  
  uint8_t* seqlib_add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q) const;
  
  uint8_t* seqlib_make_pac(const UnalignedSequenceVector& v, bool for_only) const;
  
  bwt_t *seqlib_bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr) const;
  
  /// Populate one bntann1_t entry.  (Offset in concatenated reference, etc.)
  bntann1_t* seqlib_add_to_anns(const std::string& name, const std::string& seq, bntann1_t* ann, size_t offset) const;

  friend class BWAAligner;
};

 using BWAIndexPtr = std::shared_ptr<BWAIndex>;
 
} // namespace SeqLib


