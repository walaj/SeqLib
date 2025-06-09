#pragma once

#include <cstdint>
#include <string>
#include <memory>
#include <cassert>
#include <cstdlib>
#include <unistd.h>      // for access()
#include <iostream>
#include <cstring>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>

namespace SeqLib {

// Custom deleters for HTSlib types
struct Bam1Deleter   { void operator()(bam1_t*    x) const { bam_destroy1(x); } };
struct BamHdrDeleter { void operator()(bam_hdr_t*  x) const { bam_hdr_destroy(x); } };
struct HtsFileDeleter{ void operator()(htsFile*    x) const { if(x) sam_close(x); } };
struct HtsIdxDeleter { void operator()(hts_idx_t* x) const { if(x) hts_idx_destroy(x); } };
struct HtsItrDeleter { void operator()(hts_itr_t* x) const { if(x) hts_itr_destroy(x); } };
struct BgzfDeleter   { void operator()(BGZF*      x) const { if(x) bgzf_close(x); } };

// Smart-pointer aliases
using Bam1Ptr    = std::shared_ptr<bam1_t>;
using BamHdrPtr  = std::unique_ptr<bam_hdr_t, BamHdrDeleter>;
using HtsFilePtr = std::unique_ptr<htsFile,   HtsFileDeleter>;
using HtsIdxPtr  = std::unique_ptr<hts_idx_t,  HtsIdxDeleter>;
using HtsItrPtr  = std::unique_ptr<hts_itr_t,  HtsItrDeleter>;
using BgzfPtr    = std::unique_ptr<BGZF,      BgzfDeleter>;

// Counter for BamReader progress
struct ReadCount {
    uint32_t keep{0};   ///< reads kept
    uint32_t total{0};  ///< reads seen

    // Percent of reads kept (0-100)
    int percent() const {
        return total ? static_cast<int>(100.0 * keep / total) : 0;
    }
    // Comma-formatted strings
    //std::string totalString() const { return AddCommas<uint64_t>(total); }
    //std::string keepString()  const { return AddCommas<uint64_t>(keep);  }
};

} // namespace SeqLib


// #ifndef SEQLIB_BAM_WALKER_H__
// #define SEQLIB_BAM_WALKER_H__

// #include <cassert>

// #include <stdint.h> 
// #include "SeqLib/BamRecord.h"

// // not sure what going on here...
// #ifndef INT32_MAX
// #define INT32_MAX 0x7fffffffL
// #endif

// struct idx_delete {
//   void operator()(hts_idx_t* x) { if (x) hts_idx_destroy(x); }
// };

// struct hts_itr_delete {
//   void operator()(hts_itr_t* x) { if (x) hts_itr_destroy(x); }
// };

// struct bgzf_delete {
//   void operator()(BGZF* x) { if(x) bgzf_close(x); }
// };

// struct bam_hdr_delete {
//   void operator()(bam_hdr_t* x) { if (x) bam_hdr_destroy(x); }
// };

// struct htsFile_delete { // shoudl also close cram index
//   void operator()(htsFile* x) { if (x) sam_close(x); }
// };

// // Phred score transformations
// inline int char2phred(char b) {
//   uint8_t v = b;
//   assert(v >= 33);
//   return v - 33;
// }

// // from samtools
// inline char *samfaipath(const char *fn_ref)
// {
//   char *fn_list = 0;
//   if (fn_ref == 0) return 0;
//   fn_list = (char*)calloc(strlen(fn_ref) + 5, 1);
//   strcat(strcpy(fn_list, fn_ref), ".fai");
//   if (access(fn_list, R_OK) == -1) { // fn_list is unreadable
//     std::cerr << "ERROR: Cannot read the index file for CRAM read/write" << std::endl;
//   }
//   return fn_list;
// }

// namespace SeqLib {

//   /** Small class to store a counter to measure BamReader progress.
//    * Currently only stores number of reads seen / kept. 
//    */
// struct ReadCount {

//   uint32_t keep; ///< Store total number of reads kept
//   uint32_t total; ///< Store total number of reads seen

//   ReadCount() : keep(0), total(0) {}
  
//   /** Return the percent of total reads kept
//    */
//   int percent () const {
//     int perc  = SeqLib::percentCalc<uint64_t>(keep, total); 
//     return perc;
//   }

//   /** Return the total reads visited as a comma-formatted string */
//   std::string totalString() const {
//     return SeqLib::AddCommas<uint64_t>(total);
//   }

//   /** Return the kept reads as a comma-formatted string */
//   std::string keepString() const {
//     return SeqLib::AddCommas<uint64_t>(keep);
//   }

// };

// }
// #endif
