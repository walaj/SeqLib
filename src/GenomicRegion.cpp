/*
 * By downloading the PROGRAM you agree to the following terms of use:
 *
 * BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
 * FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
 *
 * This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
 * WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
 * WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
 * NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
 *
 * 1. DEFINITIONS
 * 1.1"PROGRAM" shall mean copyright in the object code and source code known as Swap and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.or/cancer/cga/mutect on the EFFECTIVE DATE.
 *
 * 2. LICENSE
 * 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
 * LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
 * The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
 * 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
 * 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
 *
 * 3. OWNERSHIP OF INTELLECTUAL PROPERTY
 * LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
 *
 * Copyright 2012 Broad Institute, Inc.
 * Notice of attribution:  The Swap program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc. and is published at doi: 10.1038/nbt.2514.
 *
 * LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
 *
 * 4. INDEMNIFICATION
 * LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
 *
 * 5. NO REPRESENTATIONS OR WARRANTIES
 * THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
 * IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 *
 * 6. ASSIGNMENT
 * This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
 *
 * 7. MISCELLANEOUS
 * 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
 * 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
 * 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
 * 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
 * 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
 * 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
 * 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
 */

#include "SnowTools/GenomicRegion.h"

#include <cassert>
#include "SnowTools/gzstream.h" 

namespace SnowTools {

// return the width of the genomic region
int GenomicRegion::width() const {
  return pos2 - pos1 + 1;
}

// returns 0 for no overlaps, 1 for partial and 2 for complete
int GenomicRegion::getOverlap(const GenomicRegion gr) const {

  if (gr.chr != chr)
    return 0;
  
  bool gr1_in = gr.pos1 >= pos1 && gr.pos1 <= pos2;
  bool gr2_in = gr.pos2 >= pos1 && gr.pos2 <= pos2;
  bool pos1_in = pos1 >= gr.pos1 && pos1 <= gr.pos2;
  bool pos2_in = pos2 >= gr.pos1 && pos2 <= gr.pos2;

  if ( (gr1_in && gr2_in) || (pos1_in && pos2_in) )
    return 2;

  if (gr1_in || gr2_in || pos1_in || pos2_in)
    return 1;

  return 0;

}

// write genomic region to a string
std::string GenomicRegion::toString() const {
  std::stringstream out;
  //out << chrToString(chr)  << ":" << SnowUtils::AddCommas<int>(pos1) << "-" << SnowUtils::AddCommas<int>(pos2) << "(" << strand << ")"; 
  out << chrToString(chr)  << ":" << SnowTools::AddCommas<int>(pos1) << "-" << AddCommas<int>(pos2) << "(" << 
    (strand ? "+" : "-") << ")"; 
  return out.str();
}

void GenomicRegion::pad(uint32_t pad) {
  if (pad > pos1)
    pos1 = 1;
  else
    pos1 = pos1-pad;

  const uint32_t maxpos = 250000000;
  pos2 = std::min(pos2+pad, maxpos); // 2500000000 is dummy for now. should be chr end
}

bool GenomicRegion::operator<(const GenomicRegion& b) const {
  return (chr < b.chr) || (chr == b.chr && pos1 < b.pos1) || (chr==b.chr && pos1 == b.pos1 && pos2 < b.pos2);
}

bool GenomicRegion::operator==(const GenomicRegion &b) const {
  return (chr == b.chr && pos1 == b.pos1 && b.pos2 == pos2);
}

bool GenomicRegion::operator<=(const GenomicRegion &b) const {
  return (*this < b || *this == b);
}

std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr) {
  out << gr.toString();
  return out;
}

// constructor for SnowTools::GenomicRegion that takes strings. Assumes chr string is in 
// natural (1, ..., X) or (chr1, ..., chrX) format. That is, it converts to
// BamTools format with a -1 operation.
GenomicRegion::GenomicRegion(std::string t_chr, std::string t_pos1, std::string t_pos2) {

  chr = GenomicRegion::chrToNumber(t_chr);
  try {
    t_pos1 = SnowTools::scrubString(t_pos1, ",");
    t_pos2 = SnowTools::scrubString(t_pos2, ",");
    pos1 = std::stoi(t_pos1);
    pos2 = std::stoi(t_pos2);
  } catch (...) { 
    std::cerr << "stoi failed in GenomicRegion constructor. Tried: " << t_pos1 << " " << t_pos2 << std::endl;
  }
}

// constructor to take a pair of coordinates to define the genomic interval
GenomicRegion::GenomicRegion(int32_t t_chr, uint32_t t_pos1, uint32_t t_pos2, bool t_strand) {
  chr = t_chr;
  pos1 = t_pos1;
  pos2 = t_pos2;
  strand = t_strand;
}

// convert a chromosome string into a number
int GenomicRegion::chrToNumber(std::string ref) {

  // remove the chr identifier if it is there
  if (ref.find("chr") != std::string::npos)
    ref = ref.substr(3, ref.size() - 3);

  std::string ref_id = ref;
  if (ref_id == "X")
    ref_id = "23";
  else if (ref_id == "Y")
    ref_id = "24";
  else if (ref_id == "M" || ref_id == "MT")
    ref_id = "25";
  
  int out = -1;
  try {
    out = std::stoi(ref_id);
  } catch (...) {
    //cerr << "Caught error trying to convert " << ref << " to number" << endl;
  }

  //assert(out > 0);
  return (out-1); // offset by one becuase chr1 = 0 in BamAlignment coords
}

// convert a chromosome number to a string. Assumes 
// a natural ordering (1, ...), not BamTools ordering (0, ...)
std::string GenomicRegion::chrToString(int32_t ref) {
  std::string ref_id;
  if (ref == 22)
    ref_id = "X";
  else if (ref == 23)
    ref_id = "Y";
  else if (ref == 24)
    ref_id = "M";
  else
    ref_id = std::to_string(ref+1);
  assert(ref_id != "23");
  return ref_id;
}

// checks whether a GenomicRegion is empty
bool GenomicRegion::isEmpty() const {
  return chr == 0 && pos1 == 0 && pos2 == 0;
}


//
/*uint32_t GenomicRegion::posToBigPos(int refid, int pos) {
  
  if (refid < 25)
    return 0;
  
  return CHR_CLEN[refid] + pos;
  
  }*/


int GenomicRegion::distance(const GenomicRegion &gr) const {

  if (gr.chr != chr)
    return -1;
  else
    return ((pos1 > gr.pos1) ? (pos1 - gr.pos1) : (gr.pos1 - pos1));

}

void GenomicRegion::random() {
  
  uint32_t big;
  SnowTools::genRandomValue(big, SnowTools::genome_size_XY);
  
  for (size_t k = 0; k < 25; k++)
    if (big < SnowTools::CHR_CLEN[k]) {
      assert(k > 0);
      chr = --k;
      assert(big > SnowTools::CHR_CLEN[chr]);
      pos1 = big - SnowTools::CHR_CLEN[chr];
      pos2 = pos1;
      return;
    }
  std::cerr << "Value of " << big << " outside of expected range."  << std::endl;
  exit(EXIT_FAILURE);
  
}

}
