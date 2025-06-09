#pragma once

#include <string>
#include <optional>

#include "SeqLib/BamWalker.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegionCollection.h"

namespace SeqLib {

/**
 * BamReader: single-BAM reader with support for walking one or more regions.
 * Provides a pull-style Next() API returning std::optional<BamRecord>.
 */
class BamReader {
public:
  BamReader() = default;
  ~BamReader() = default;

  /**
   * Open a BAM/SAM/CRAM file for reading.
   * @param path filesystem path or "-" for stdin
   * @return true on success
   */
  bool Open(const std::string& path);

  /** Close the file and clear all state. */
  void Close();

  /** Reset back to the start of the file or current regions. */
  void Reset();

  /**
   * Restrict iteration to a single region. Clears any prior regions.
   * @return true if iterator is successfully positioned
   */
  bool SetRegion(const GenomicRegion& region);

  /**
   * Restrict iteration to multiple regions in order. Clears prior regions.
   * @return true if first region iterator is valid
   */
  bool SetRegions(const GRC& regions);

  /**
   * Pull the next record. Returns std::nullopt at EOF or end of regions.
   */
  std::optional<BamRecord> Next();

  /** Access the BAM header. */
  const BamHeader& Header() const;

  /** True if a file is currently open. */
  bool IsOpen() const { return static_cast<bool>(fp_); }

  /** Set the path to the reference needed to unpack a CRAM file */
  void SetCramReference(const std::string& cram);

  /** Print the regions this BAM will iterate over */
  std::string PrintRegions() const;

  friend std::ostream& operator<<(std::ostream& out, const BamReader& b);
  
private:
  HtsFilePtr fp_;         ///< HTSlib file handle
  HtsIdxPtr  idx_;        ///< Index for random access
  HtsItrPtr  itr_;        ///< Iterator for current region
  BamHeader  hdr_;        ///< SeqLib Bam Header
  
  GRC        regions_;    ///< List of regions to walk
  size_t     region_idx_ = 0; ///< Current index in regions_

  std::string path_;
  std::string cram_reference_; 
};

} // namespace SeqLib
