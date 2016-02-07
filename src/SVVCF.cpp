#include "SnowTools/SVVCF.h"

#include "SnowTools/SnowUtils.h"
#include "SnowTools/gzstream.h"
#include <stdexcept>
#include <boost/regex.hpp>

namespace SnowTools {


  SVVCF::SVVCF(const std::string& vcf_file, bam_hdr_t *h) {

    // make sure it exists
    if (!read_access_test) 
      std::runtime_error("SnowTools::SVVCF Could not open file: " + vcf_file);

    igzstream iss(vcf_file.c_str());
    if (!iss || vcf_file.length() == 0) 
      std::runtime_error("SnowTools::SVVCF Could not open file: " + vcf_file);
    
    std::string line;
    while (std::getline(iss, line, '\n')) {
      if (line.length() > 0) {

	if (line.at(0) == '#') // ignore comments
	  continue;

	std::istringstream iss_this(line);
	int count = 0;
	std::string val, chr, pos, id, ref, alt, info;
	
	while (std::getline(iss_this, val, '\t')) {
	  switch (count) {
	  case 0 : chr = val;
	  case 1 : pos = val;
	  case 2 : id = val;
	  case 3 : ref = val;
	  case 4 : alt = val;
	  case 7: info = val;
	  }
	  ++count;
	}

	if (count < 3) {
	  std::cerr << "Didn't parse VCF line properly: " << line << std::endl;
	  return;
	}
	
	// parse out the id
	boost::regex reg_id("([0-9]+):.*"); 
	boost::cmatch match_id;
	if (boost::regex_match(id.c_str(), match_id, reg_id)) {
	  try { 
	    id = match_id[1].str();
	  } catch (...) { 
	    std::cerr << "Can't parse ID for " << id << std::endl;
	  }
	}

	// already added?
	if (m_hash.count(id))
	  continue;

	// parse out insertion
	std::string ins;
	boost::regex reg_ins(".*INSERTION=([A-Z]+).*");
	boost::cmatch match_ins;
	if (boost::regex_match(info.c_str(), match_ins, reg_ins)) {
	  ins = match_ins[1].str();
	}

	// parse out homology
	std::string hom;
	boost::regex reg_hom(".*HOMSEQ=([A-Z]+).*");
	boost::cmatch match_hom;
	if (boost::regex_match(info.c_str(), match_hom, reg_hom)) {
	  hom = match_hom[1].str();
	}

	std::string chr2, pos2;

	// parse out the alt location
	boost::regex reg(".*?(\\[|\\])(.*?):([0-9]+).*"); 
	boost::cmatch match;
	if (boost::regex_match(alt.c_str(), match, reg)) {
	  try { 
	    chr2 = match[2].str();
	    pos2 = match[3].str();
	  } catch (...) { 
	    std::cerr << "Can't parse ALT line for " << alt << std::endl;
	  }
	}
	
	// get strandedness (Marcin convention, + is direction of DNA, away from break)
        char strand2 = (alt.back() == '[' || alt.at(0) == '[') ? '+' : '-';
	char strand1 = (alt.at(0) == '[' || alt.at(0) == ']') ? '+' : '-';

	// construct the GenomicRegion
	GenomicRegion gr1(chr, pos, pos, h);
	GenomicRegion gr2(chr2, pos2, pos2, h);
	gr1.strand = strand1;
	gr2.strand = strand2;
	if (gr1.valid() && gr2.valid()) {
	  m_hash.insert(std::pair<std::string, int>(id, m_side1.size()));
	  m_side1.add(gr1);
	  m_side2.add(gr2);
	  m_insertions.push_back(ins);
	  m_homology.push_back(hom);
	  m_name.push_back(id);
	  assert(m_hash.size() == m_side1.size());
	  assert(m_side1.size() == m_side2.size());
	  assert(m_insertions.size() == m_side2.size());
	} else {
	  std::cerr << "SVVCF - Couldn't parse correctly " << line << std::endl;
	  std::cerr << " _____________ " << gr1 << " -- " << gr2 << std::endl;
	}
	
      }
    } // end while
    
    
  }

  size_t SVVCF::size() const {
    return m_hash.size();
  }

  std::string SVVCF::getSurroundingContigFromReference(size_t i, RefGenome& ref, const bam_hdr_t *h, BamReadVector& brv, int width) const {
    
    if (i >= m_side1.size())
      throw std::invalid_argument("SVVCF::getSurroundingContigFromReference - Query is out of bounds");

    const GenomicRegion * gr1 = &(m_side1[i]);
    const GenomicRegion * gr2 = &(m_side2[i]);

    int hom_len = m_homology[i].length();

    if (m_name[i] == "13929")
      std::cerr << m_homology[i] << std::endl;

    int p1 = gr1->strand == '+' ? gr1->pos1 + hom_len : gr1->pos1 - width;
    int p2 = gr1->strand == '+' ? gr1->pos1 + width : gr1->pos1;
    int d1 = gr2->strand == '+' ? gr2->pos1 : gr2->pos1 - width;
    int d2 = gr2->strand == '+' ? gr2->pos1 + width : gr2->pos1;

    std::string seq1 = ref.queryRegion(gr1->ChrName(h), p1, p2);
    std::string seq2 = ref.queryRegion(gr2->ChrName(h), d1, d2);
    
    // make the bam reads
    std::string name = "id_" + m_name[i]; //std::to_string(i);

    // when making contig, need to reverse one
    // set convention that (1) is always true to gr1 orientation 
    // (if gr1 is (-) then put on rev strand)
    // ++ : pos1 rev2    g : +-     // g is orientation for bam flag
    // +- : pos1 pos2    g : ++
    // -+ : pos1 pos2    g : ++
    // -- : rev1 pos2    g : -+
    GenomicRegion g1 = *gr1;
    GenomicRegion g2 = *gr2;

    bool inv = gr1->strand == gr2->strand;
    if (inv && gr1->strand == '-') // --
      { rcomplement(seq1); g2.strand = '+'; }
    else if (inv && gr2->strand == '+') // ++
      { rcomplement(seq2); g2.strand = '-'; }
    else if (!inv && gr2->strand == '-') // +-
      g2.strand = '+';
    else if (!inv && gr2->strand == '+') // -+
      g1.strand = '+';

    std::string seq = seq1 + m_insertions[i] + seq2;
    int is = m_insertions[i].length();

    // make the cigar strings for the two (soft clip the other one)
    Cigar cig1, cig2;
    uint32_t c1_1, c1_2, c2_1, c2_2;
    c1_1 = seq1.length() << BAM_CIGAR_SHIFT | BAM_CMATCH;
    c1_2 = (seq2.length() + is) << BAM_CIGAR_SHIFT | BAM_CSOFT_CLIP;
    c2_1 = (seq1.length() + is) << BAM_CIGAR_SHIFT | BAM_CSOFT_CLIP;
    c2_2 = seq2.length() << BAM_CIGAR_SHIFT | BAM_CMATCH;
    
    cig1.push_back(CigarField(c1_1));
    cig1.push_back(CigarField(c1_2));
    cig2.push_back(CigarField(c2_1));
    cig2.push_back(CigarField(c2_2));

    if (inv && gr1->strand == '-')  {// flip it back to positive (--)
      rcomplement(seq);
      std::reverse(cig1.begin(), cig1.end());
    }
    BamRead b1(name, seq, &g1, cig1);
    if (inv && gr1->strand == '-')  // put it back (--)
      rcomplement(seq);      

    if (inv && gr2->strand == '+') { // if ++, then (2) was flipped to (-). Put back to 
      rcomplement(seq);
      std::reverse(cig2.begin(), cig2.end());
    }
    BamRead b2(name, seq, &g2, cig2);
    if (inv && gr2->strand == '+') // put it back (++)
      rcomplement(seq);

    b1.AddZTag("MC", std::to_string(gr1->chr+1));
    b2.AddZTag("MC", std::to_string(gr2->chr+1));
    
    brv.push_back(b1);
    brv.push_back(b2);

    return seq;
    
  }
  
}
