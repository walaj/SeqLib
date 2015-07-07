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
 * 1.1"PROGRAM" shall mean copyright in the object code and source code known as MuTect and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.or/cancer/cga/mutect on the EFFECTIVE DATE.
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
 * Notice of attribution:  The MuTect program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc. and is published at doi: 10.1038/nbt.2514.
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
v * 7. MISCELLANEOUS
 * 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
 * 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
 * 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
 * 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
 * 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
 * 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
 * 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
 */

#include "SnowTools/Histogram.h"
#include "SnowTools/SnowUtils.h"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>

#define BINARY_SEARCH 1

#define DEBUG_HISTOGRAM

namespace SnowTools {

Histogram::Histogram(const int32_t& start, const int32_t& end, const uint32_t& width)
{
  
  assert(end >= start);

  Bin bin;
  bin.bounds.first = start;

  int32_t next_end = start + width - 1; // -1 because width=1 is bound.first = bounds.second

  while (next_end < end) 
    {
      // finish this bin
      bin.bounds.second = next_end;
      m_bins.push_back(bin);
      m_ind.push_back(bin.bounds.first); // make the index for the lower bound

      // start a new one
      bin.bounds.first = next_end+1;
      next_end += width;
    }
  
  // finish the last bin
  bin.bounds.second = end;
  m_bins.push_back(bin);
  m_ind.push_back(bin.bounds.first);

  // add a final bin
  bin.bounds.first = end+1;
  bin.bounds.second = 250000000;
  m_bins.push_back(bin);
  m_ind.push_back(bin.bounds.first);

}

void Histogram::toCSV(std::ofstream &fs) {

  for (auto& i : m_bins) 
    fs << i << std::endl;

}
void Histogram::removeElem(const int32_t& elem) {
  --m_bins[retrieveBinID(elem)];
}

void Histogram::addElem(const int32_t& elem) {
  ++m_bins[retrieveBinID(elem)];
}

std::string Histogram::toFileString() const {
  std::stringstream ss;
  for (auto& i : m_bins)
    if (i.m_count)
      ss << i.bounds.first << "_" << i.bounds.second << "_" << i.m_count << ",";
  return(cutLastChar(ss.str())); // trim off last comma
  
}

size_t Histogram::retrieveBinID(const int32_t& elem) const {

  if (elem < m_bins[0].bounds.first) 
    {
#ifdef DEBUG_HISTOGRAM
      std::cerr << "removeElem: elem of value " <<  elem << " is below min bin " << m_bins[0] << std::endl;
      exit(1);
#endif
      return 0;
    }

  if (elem > m_bins.back().bounds.second) 
    {
#ifdef DEBUG_HISTOGRAM
      std::cerr << "removeElem: elem of value " <<  elem << " is above max bin " << m_bins.back() << std::endl;
      exit(1);
#endif
      return m_bins.size();
    }


  if (m_bins[0].contains(elem)) 
    return 0;
  if (m_bins.back().contains(elem)) 
    return m_bins.size();

#ifdef BINARY_SEARCH
  // binary search
  std::vector<int32_t>::const_iterator it = std::upper_bound(m_ind.begin(), m_ind.end(), elem);
  size_t i = it - m_ind.begin()-1;
  assert(i < m_ind.size());
  return i;
#else
  for (size_t i = 0; i < m_bins.size(); i++) {
    if (m_bins[i].contains(elem)) {
      return i;
    }
  }
#endif
  std::cerr << "bin not found for element " << elem << std::endl;
  return 0;
}

void Histogram::initialSpans(size_t num_bins, std::vector<S>* pspanv, size_t min_bin_width) {

  // ensure that they spans are sorted
  std::sort(pspanv->begin(), pspanv->end());

  // fill the histogram bins with matrix pairs (pre-sorted by distance)
  Bin bin; 

  // get number of inter-chr events
  size_t intra = 0;
  for (auto& i : *pspanv)
    if (i != INTERCHR)
      intra++;

  size_t bin_cut = 0;
  try {
    bin_cut = floor((double)intra / (double)num_bins);
    if (bin_cut == 0)
      throw 1;
  } catch(...) {
    std::cerr << "Error in determining bin cut. Not enought events or too many bins?" << std::endl;
    std::cerr << "Events: " << pspanv->size() << " Num Bins " << num_bins << " quantile count (hist height) " << bin_cut << std::endl;
  }

  std::cout << "...Events per bin: " << bin_cut << " num bins " << num_bins << std::endl;

  S last_span = 0;
  size_t tcount = 0; // count events put into bins

  // iterate over spans
  for (auto& span : *pspanv) {
    if (span != INTERCHR) {
      
      ++tcount;
      
      // moved into a new bin? (or done?)
      if (bin.getCount() > bin_cut && span != last_span && (last_span - bin.bounds.first) >= min_bin_width) { 

	// finalize, save old bin
	bin.bounds.second = last_span;
	m_bins.push_back(bin);
	
	// new bin
	bin.bounds.first = last_span+1;
	bin.m_count = 0;
	
      }
      ++bin;
      if (bin.getCount() >= bin_cut) {
	last_span = span;
      }
      
      //update the size of current bin
      bin.bounds.second = span;
    }
  }
  // add the last bin
  bin.bounds.second = INTERCHR-1; 
  m_bins.push_back(bin);

  // add a bin for interchr events
  bin.bounds = {INTERCHR, INTERCHR};
  bin.m_count = pspanv->size() - intra;
  m_bins.push_back(bin);

  // make the indices of lower bound
  for (auto& i : m_bins)
    m_ind.push_back(i.bounds.first);

  if (m_bins.size() != (num_bins+1)) {
    std::cout << " bin cut " << bin_cut << std::endl;
    std::cout << " num bins " << num_bins << " bins.size() " << m_bins.size() << std::endl;
    //assert(bins.size() == (num_bins+1));
  }

}

bool Bin::operator < (const Bin& b) const {
  return (bounds.first < b.bounds.first || (bounds.first==b.bounds.first && bounds.second < b.bounds.second));

}

bool Bin::contains(const int32_t& elem) const {

  return (elem >= bounds.first && elem <= bounds.second); 


}

Bin& Bin::operator++()
{
  ++m_count;
  return *this;
}


Bin& Bin::operator--() {
  assert(m_count > 0); 
  --m_count;
  return *this;
}

}
