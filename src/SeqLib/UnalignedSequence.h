#ifndef SEQLIB_UNALIGNED_SEQ_H__
#define SEQLIB_UNALIGNED_SEQ_H__

#include <cstring>

namespace SeqLib {

  /** Structure to hold unaligned sequence (name and bases)
   */
  struct UnalignedSequence {
    std::string Name; ///< Name of the contig
    std::string Seq; ///< Sequence of the contig (upper-case ACTGN)
    std::string Qual; ///< Quality scores
  };
  typedef std::vector<UnalignedSequence> UnalignedSequenceVector; ///< A collection of contigs that can be indexed

}

#endif
