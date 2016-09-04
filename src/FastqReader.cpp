#include "SeqLib/FastqReader.h"

namespace SeqLib {

  FastqReader::FastqReader(const std::string& file) : m_iss(NULL) {

  m_file = file;

  // check if file exists
  struct stat buffer;   
  if (stat (file.c_str(), &buffer) != 0) {
    std::cerr << "FastqReader: Failed to read non-existant file " << m_file << std::endl;
    return;
  }

  // open the file
  m_iss = new std::ifstream(file.c_str());

  if (!m_iss) {
    std::cerr << "FastqReader: Failed to read " << m_file << std::endl;
    return;
  }

}

bool FastqReader::GetNextSequence(std::string& qn, std::string& seq) {

  if (!m_iss)
    return false;

  std::string line;
  size_t c = 0;

  // loop through the lines until run out
  while(std::getline(*m_iss, line, '\n')) {
    ++c;

    // qname
    if (c == 1) {
      if (line.length() < 2) {
	std::cerr << "FastqReader: Empty qname on line " << line << std::endl;
	return false;
      }

      // get the qname without the leading @ sign
      qn = line.substr(1,line.length() - 1);
    } else if (c == 2) {
      seq = line;
      return true; // fasta reader
    } else if (c == 4) { // we're at the end
      return true;
    }
    
  }

  // we didn't get into the loop, so we are out of reads
  return false;

}

}
