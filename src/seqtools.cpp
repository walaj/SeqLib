#include <getopt.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>

#include "SeqLib/BFC.h"
#include "SeqLib/FastqReader.h"

#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

static const char *SEQTOOLS_USAGE_MESSAGE =
"Program: seqtools \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: seqtools <command> [options]\n\n"
"Commands:\n"
"           bfc       \n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

void runbfc(int argc, char** argv);
void parseBfcOptions(int argc, char** argv);

namespace bfcopt {

  static std::string input;
}

static const char* shortopts = "hi:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input",   required_argument, NULL, 'i' },
  { NULL, 0, NULL, 0 }
};

int main(int argc, char** argv) {

  if (argc <= 1) {
    std::cerr << SEQTOOLS_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << SEQTOOLS_USAGE_MESSAGE;
      return 0;
    } else if (command == "bfc") {
      runbfc(argc -1, argv + 1);
    } else {
      std::cerr << SEQTOOLS_USAGE_MESSAGE;
      return 0;
    }
  } 
  
  return 0;

}

void runbfc(int argc, char** argv) {

  parseBfcOptions(argc, argv);

  // read in a fasta file
  SeqLib::FastqReader f(bfcopt::input);

  SeqLib::BFC b;

  std::string qn, seq;
  while (f.GetNextSequence(qn, seq)) {
    std::string e;
    assert(b.AddSequence(seq.c_str(), e.c_str(), qn.c_str()));
  }


  assert(b.Train());
  assert(b.ErrorCorrect());

  SeqLib::UnalignedSequenceVector u;
  b.GetSequences(u);
  std::cerr << " U " << u.size() << std::endl;
  std::cerr << " KCOV " << b.GetKCov() << std::endl;
  std::cerr << " KMER " << b.GetKMer() << std::endl;

  for (SeqLib::UnalignedSequenceVector::const_iterator i = u.begin();
       i != u.end(); ++i) {
    std::cout << ">" << i->Name << std::endl << i->Seq << std::endl;
  }

  
  
}

// parse the command line options
void parseBfcOptions(int argc, char** argv) {

  bool die = false;
  bool help = false;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'i': arg >> bfcopt::input; break;
    default: die= true; 
    }
  }

  if (die || help) {
      std::cerr << "\n" << SEQTOOLS_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}

