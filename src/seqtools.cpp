#include <getopt.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>

#include "SeqLib/BFC.h"
#include "SeqLib/FastqReader.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BWAWrapper.h"

#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

static const char *SEQTOOLS_USAGE_MESSAGE =
"Program: seqtools \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: seqtools <command> [options]\n\n"
"Commands:\n"
"           bfc       \n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

static const char *BFC_USAGE_MESSAGE =
"Program: seqtools bfc \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: seqtools bfc [options]\n\n"
"Commands:\n"
"  --verbose,   -v        Set verbose output\n"
"  --fasta,     -f        Output stream should be a FASTA (no realignment)\n"
"  --bam,       -b        Output stream should be a BAM (not SAM)\n"
"  --cram,      -C        Output stream should be a CRAM (not SAM)\n"
"  --infasta,   -F <file> Input a FASTA insted of BAM/SAM/CRAM stream\n"
"  --reference, -G <file> Reference genome if using BWA-MEM realignment\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

void runbfc(int argc, char** argv);
void parseBfcOptions(int argc, char** argv);

namespace opt {

  static bool verbose = false;
  static char mode = 's';
  static std::string input;
  static std::string reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  static std::string fasta; // input is a fasta
}

static const char* shortopts = "hbfvCG:F:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "bam",                     no_argument, NULL, 'b' },
  { "cram",                    no_argument, NULL, 'C' },
  { "fasta",                   no_argument, NULL, 'f' },
  { "infasta",                 required_argument, NULL, 'F' },
  { "reference",               required_argument, NULL, 'G' },
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

  SeqLib::BFC b;

  // is this a fasta file
  
  if (!opt::fasta.empty()) {
    // read in a fasta file
    SeqLib::FastqReader f(opt::fasta);
    
    std::string qn, seq;
    while (f.GetNextSequence(qn, seq)) {
      std::string e;
      assert(b.AddSequence(seq.c_str(), e.c_str(), qn.c_str()));
    }
  } else { //if (opt::mode == 'b' || opt::mode == 's' || opt::mode == 'C') {
    SeqLib::BamReader br;
    br.Open(opt::input == "-" ? "-" : opt::input);
    SeqLib::BamRecord rec;
    while(br.GetNextRecord(rec)) {
      b.AddSequence(rec); //rec.Sequence().c_str(), rec.Qualities().c_str(), rec.Qname().c_str());
    }
  } 

  if (!b.Train()) {
    std::cerr << "Training failed on " << b.NumSequences() << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!b.ErrorCorrect()) {
    std::cerr << "Correction failed on " << b.NumSequences() << std::endl;
    exit(EXIT_FAILURE);
  }

  SeqLib::UnalignedSequenceVector u;
  b.GetSequences(u);
  if (opt::verbose)
    std::cerr << "nseqs: " << u.size() 
	      << " kcov: " << b.GetKCov() 
	      << " kmer: " << b.GetKMer() << std::endl;
  
  
  if (opt::mode == 'f') {
    for (SeqLib::UnalignedSequenceVector::const_iterator i = u.begin();
	 i != u.end(); ++i) {
      std::cout << ">" << i->Name << std::endl << i->Seq << std::endl;
    }
    return;
  } 

  SeqLib::BamWriter bw;
  if (opt::mode == 'b')
    bw = SeqLib::BamWriter(SeqLib::BAM);
  else if (opt::mode == 's')
    bw = SeqLib::BamWriter(SeqLib::SAM);
  else if (opt::mode == 'C') {
    bw = SeqLib::BamWriter(SeqLib::CRAM);
    bw.SetCramReference(opt::reference);
  }
  else {
    std::cerr << "Unrecognized output stream mode " << opt::mode << std::endl;
    exit(EXIT_FAILURE);
  }
  
  bw.Open("-");
  
  SeqLib::BWAWrapper bwa;
  if (!bwa.LoadIndex(opt::reference)) {
    std::cerr << "Failed to load index for BWA-MEM from: " << opt::reference << std::endl;
    exit(EXIT_FAILURE);
  }
  
  bw.SetHeader(bwa.HeaderFromIndex());
  bw.WriteHeader();
  
  // run through and read
  for (SeqLib::UnalignedSequenceVector::const_iterator i = u.begin(); i != u.end(); ++i) {
    SeqLib::BamRecordVector brv;
    bool hardclip = false;
    double frac = 0.9;
    int max_secondary = 10;
    bwa.AlignSequence(i->Seq, i->Name, brv, false, frac, 10);
    for (SeqLib::BamRecordVector::iterator r = brv.begin();
	 r != brv.end(); ++r) {
      if (!i->Qual.empty())
	r->SetQualities(i->Qual, 33);
      bw.WriteRecord(*r);
    }
  }
  
}

// parse the command line options
void parseBfcOptions(int argc, char** argv) {

  bool die = false;
  bool help = false;

  // get the first argument as input
  if (argc > 1)
    opt::input = std::string(argv[1]);
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'f': opt::mode = 'f'; break;
    case 'F': arg >> opt::fasta; break;
    case 'b': opt::mode = 'b'; break;
    case 'C': opt::mode = 'C'; break;
    case 'G': arg >> opt::reference; break;
    default: die= true; 
    }
  }

  if (die || help || (opt::input.empty() && opt::fasta.empty())) {
      std::cerr << "\n" << BFC_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}

