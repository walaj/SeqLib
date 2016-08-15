#ifndef SNOWTOOLS_BAM_HEADER_H__
#define SNOWTOOLS_BAM_HEADER_H__

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"

#include <string>
#include <memory>
#include <unordered_map>

namespace SeqKit {
  
  /** Store a header to a BAM file 
   *
   * Stores a BAM header, which also acts as a dictionary of 
   * reference sequences, with names and lengths.
   */
  class BamHeader {

  public:

    /** Initializes a new empty BamHeader with no data
     * 
     * @note No memory is allocated here
     */
    BamHeader() {};
    
    /** Initialize a BamHeader from a string containing
     * a BAM header in human-readable form (e.g. @PG ... )
     * @param Text of a BAM header, with newlines separating lines
     */
    BamHeader(const std::string& hdr);

    /** Create a new BamHeader from a raw HTSlib header.
     * 
     * @note This will make a copy of the input header
     */
    BamHeader(const bam_hdr_t * hdr);
    
    /** Return the number of sequences store in this dictionary
     * Returns 0 if header is unitialized.
     */
    int NumSequences() const;

    /** Convert a numeric sequence ID to a name
     * 
     * @exception Throws an out_of_range if ID is >= then number of 
     * targets in dictionary, or if header is uninitialized..
     * @exception Throws an invalid_argument if ID is < 0;
     */
    std::string IDtoName(int id) const;

    /** Check if the header has been initialized
     */
    bool isEmpty() const { return h.get() == 0; }

    /** Return the raw bam_hdr_t */
    const bam_hdr_t* get() const { return h.get(); }

    /** Get the numeric ID associated with a sequence name.
     * @param name Name of the sequence
     * @return ID of named sequence, or -1 if not in dictionary
     */
    int Name2ID(const std::string& name) const;

    void WriteToStdout() const;
  private:

    // adapted from sam.c - bam_nam2id
    int bam_name2id_2(const bam_hdr_t *h, const char *ref) const;

    std::shared_ptr<bam_hdr_t> h;

    // make the name 2 id map (to be used by Name2ID)
    // replaces part of bam_name2id that makes the hash table
    int ConstructName2IDTable();

    // hash table for name to id
    std::shared_ptr<std::unordered_map<std::string, int>> n2i;

    // adapted from sam_hdr_read
    bam_hdr_t* sam_hdr_read2(const std::string& hdr) const;

  };
  
}


#endif
