#include "SeqLib/FastqReader.h"

#include <cctype>
#include <algorithm>

namespace SeqLib {

  bool FastqReader::Open(const std::string& f) {

    m_file = f;

    // check if file exists
    struct stat buffer;   
    if (stat (m_file.c_str(), &buffer) != 0) {
      std::cerr << "FastqReader: Failed to read non-existant file " << m_file << std::endl;
      return false;
    }
    
    // open the file
    m_iss = new std::ifstream(m_file.c_str());
    
    if (!m_iss) {
      std::cerr << "FastqReader: Failed to read " << m_file << std::endl;
      return false;
    }

    return true;
    
  }

  FastqReader::FastqReader(const std::string& file) : m_file(file), m_type('*'), m_iss(NULL) {
    Open(m_file);
  }

bool FastqReader::GetNextSequence(UnalignedSequence& s) {

  if (!m_iss)
    return false;

  std::string line;
  size_t c = 0;

  // loop through the lines until run out
  while(std::getline(*m_iss, line, '\n')) {
    ++c;

    // stop on empty line
    if (line.empty())
      return false;

    // determine if fastq or fasta
    if (m_type == 'k') { 
      if (line.at(0) == '@')
	m_type = 'q';
      else if (line.at(0) == '>')
	m_type = 'a';
      else {
	std::cerr << "Cannot determine fasta vs fasta from line " << line << std::endl;
	return false;
      }
    }
    
    // get the qname without the leading @ sign or >
    switch (c) {
    case 1: s.Name = line.substr(1,line.length() - 1); break;
    case 2: 
      s.Seq = line; 
      transform(s.Seq.begin(), s.Seq.end(), s.Seq.begin(), ::toupper); 
      if (m_type == 'a')
	return true;
      break;
    case 3: s.Strand == line.at(0); break;
    case 4: s.Qual = line; return true;
    }
    
  }

  // we didn't get into the loop, so we are out of reads
  return false;

}

}
