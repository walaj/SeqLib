#ifndef SEQLIB_FERMI_H__
#define SEQLIB_FERMI_H__

#include <string>
#include <cstdlib>
#include <iostream>

#include "SeqLib/BamRecord.h"

extern "C" 
{
#include "fermi-lite/htab.h"
#include "fermi-lite/fml.h"
#include "fermi-lite/bfc.h"
}

namespace SeqLib {

  /** Sequence assembly using FermiKit from Heng Li
   */
  class FermiAssembler {

  public:

    /** Create an empty FermiAssembler with default parameters */
    FermiAssembler ();

    /** Destroy by clearing all stored reads from memory */
    ~FermiAssembler();

    /** Provide a set of reads to be assembled 
     * @param Reads with or without quality scores
     * @note This will copy the reads and quality scores
     * into this object. Deallocation is automatic with object
     * destruction, or with ClearReads.
     */ 
    void AddReads(const BamRecordVector& brv);

    /** Clear all of the sequences and deallocate memory.
     * This is not required, as it will be done on object destruction
     */
    void ClearReads();

    /** Clear all of the contigs and deallocate memory.
     * This is not required, as it will be done on object destruction
     */
    void ClearContigs();

    /** Peform Bloom filter error correction of the reads
     * in place. */
    void CorrectReads();

    /** Peform Bloom filter error correction of the reads
     * in place. Also remove unique reads.
     */
    void CorrectAndFilterReads();
    
    /** Return the sequences in this object, which may have 
     * been error-corrected
     */
    UnalignedSequenceVector GetSequences() const;

    /** Perform the string graph assembly.
     * This will product the string graph,
     * and travserse the graph to emit contigs
     */
    void PerformAssembly();

    /** Return the assembled contigs. 
     */
    std::vector<std::string> GetContigs() const;

    /** Perform assembly, without error correction */
    void DirectAssemble(float kcov);

    /** Set the minimum overlap between reads during string graph construction */
    void SetMinOverlap(uint32_t m) { opt.min_asm_ovlp = m; }

    /** Return the minimum overlap parameter for this assembler */
    uint32_t GetMinOverlap() const { return opt.min_asm_ovlp; }

    /** Add a set of unaligned sequences to stage for assembly */
    void AddReads(const UnalignedSequenceVector& v);

  private:

    // reads to assemble
    fseq1_t *m_seqs;
    
    std::vector<std::string> m_names;

    // number of reads
    size_t n_seqs;

    // number of contigs
    int n_utg;
  
    // options
    fml_opt_t opt;

    // the unitigs
    fml_utg_t *m_utgs;

  };
  

}

#endif
