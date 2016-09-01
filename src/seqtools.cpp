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
  //"  --input,     -i        Input FASTA, BAM, CRAM, SAM. If not specified, reads from stdin\n"
  //"  --imode,     -m        Input mode. f: FASTA  b: BAM/CRAM/SAM  <none>: stdin (sam/bam stream)\n"
  //"  --omode,     -w        Output stream mode. f: FASTA  b: BAM  s: SAM   <none>: stdin (sam/bam stream)\n"
"  --fasta,     -f        Output stream is a fasta (no realignment)"
"  --bam,       -b,       Output stream is BAM (not SAM)"
"  --reference, -G        Reference genome if using BWA-MEM realignment\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

void runbfc(int argc, char** argv);
void parseBfcOptions(int argc, char** argv);

namespace opt {

  static std::string input;
  static std::string reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  static char inputmode = 'f';
  static char outputmode = 's';
}

static const char* shortopts = "hi:m:w:G:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input",                   required_argument, NULL, 'i' },
  { "imode",                   required_argument, NULL, 'm' },
  { "omode",                   required_argument, NULL, 'w' },
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

  if (opt::inputmode == 'f') {
    // read in a fasta file
    SeqLib::FastqReader f(opt::input);
    
    std::string qn, seq;
    while (f.GetNextSequence(qn, seq)) {
      std::string e;
      assert(b.AddSequence(seq.c_str(), e.c_str(), qn.c_str()));
    }
  } else if (opt::inputmode == '-' || opt::inputmode == 'b') {
    SeqLib::BamReader br;
    br.Open(opt::input.empty() ? "-" : opt::input);
    SeqLib::BamRecord rec;
    while(br.GetNextRecord(rec)) {
      b.AddSequence(rec.Sequence().c_str(), rec.Qualities().c_str(), rec.Qname().c_str());
    }
  } else {
    std::cerr << "Input mode: " << opt::inputmode << " not recognized " << std::endl;
    exit(EXIT_FAILURE);
  }

  assert(b.Train());
  assert(b.ErrorCorrect());

  SeqLib::UnalignedSequenceVector u;
  b.GetSequences(u);
  std::cerr << "nseqs: " << u.size() 
	    << " kcov: " << b.GetKCov() 
	    << " kmer " << b.GetKMer() << std::endl;

  
  if (opt::outputmode == 'f') {
    for (SeqLib::UnalignedSequenceVector::const_iterator i = u.begin();
	 i != u.end(); ++i) {
      std::cout << ">" << i->Name << std::endl << i->Seq << std::endl;
    }
    return;
  } 

  SeqLib::BamWriter bw;
  if (opt::outputmode == 'b')
    bw = SeqLib::BamWriter(SeqLib::BAM);
  else if (opt::outputmode == 's')
    bw = SeqLib::BamWriter(SeqLib::SAM);
  else {
    std::cerr << "Unrecognized output stream mode " << opt::outputmode << std::endl;
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
      r->SetQualities(i->Qual);
      bw.WriteRecord(*r);
    }
  }
  
}

// parse the command line options
void parseBfcOptions(int argc, char** argv) {

  bool die = false;
  bool help = false;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'i': arg >> opt::input; break;
    case 'm': arg >> opt::inputmode; break;
    case 'w': arg >> opt::outputmode; break;
    case 'G': arg >> opt::reference; break;
    default: die= true; 
    }
  }

  if (die || help) {
      std::cerr << "\n" << BFC_USAGE_MESSAGE;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}

