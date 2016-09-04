#ifndef SEQLIB_FASTQ_READER_H__
#define SEQLIB_FASTQ_READER_H__

#include <string>

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>

namespace SeqLib{

class FastqReader {

 public:

  FastqReader(const std::string& file);

  bool GetNextSequence(std::string& qn, std::string& seq);

  ~FastqReader() { if (m_iss) delete m_iss; }

 private:
  
  std::string m_file;
  std::ifstream * m_iss;

};

}

#endif
