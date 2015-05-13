#include "VFilterTool.h"

#include <string>
#include <getopt.h>
#include <iostream>

#include "SnowTools/gzstream.h"
#include "SnowTools/SnowUtils.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/GenomicRegion.h"

#include "VariantBamWalker.h"

using SnowTools::GenomicRegion;
using SnowTools::GenomicRegionCollection;
using SnowTools::GRC;

static const char *VARIANT_BAM_USAGE_MESSAGE =
"Usage: vfilter -i <input.bam> -o <output.bam> [OPTIONS] \n\n"
"  Description: Process a BAM file for use with rearrangement variant callers by removing proper pairs and bad regions\n"
"\n"
" General options\n"
"  -v, --verbose                        Display more verbose output\n"
"  -h, --help                           Display this help and exit\n"
"  -i, --input-bam                      BAM file to filter\n"
"  -o, --output-bam                     Output BAM file\n"
"  -g, --region                         Regions (e.g. myvcf.vcf or WG for whole genome) or newline seperated subsequence file.  Applied in same order as -r for multiple\n"
"  -r, --rules                          A script for the rules. If specified multiple times, will be applied in same order as -g\n"
"  -f, --rules-script                   A file (script) for the rules and regions/sequences. If specified, ignores -r and -g flags.\n"
"  -k, --regions-file                   BED file of regions to process. Default: Whole genome\n"
" Optional Input\n"
"\n";

namespace opt {

  static std::string bam;
  static std::string out;
  static size_t verbose = 1;
  static std::string rules = "";
  static std::string proc_regions = "";
  static std::string refgenome = SnowTools::REFHG19;
}

static const char* shortopts = "hv:qji:o:r:c:k:g:c:f:T:";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, 'h' },
  { "verbose",                    required_argument, NULL, 'v' },
  { "input-bam",                  required_argument, NULL, 'i' },
  { "output-bam",                 required_argument, NULL, 'o' },
  { "qc-only",                    no_argument, NULL, 'q' },
  { "rules",                      required_argument, NULL, 'r' },
  { "reference",                  required_argument, NULL, 'T' },
  { "rules-script",               required_argument, NULL, 'f' },
  { "region",                     required_argument, NULL, 'g' },
  { "regions-file",          required_argument, NULL, 'k' },
  { NULL, 0, NULL, 0 }
};

// forward declare
void parseVarOptions(int argc, char** argv);

void runVFilter(int argc, char** argv) {

  // parse the command line
  parseVarOptions(argc, argv);

  if (opt::verbose > 0) {
    std::cout << "Input BAM:  " << opt::bam << std::endl;
    std::cout << "Output BAM: " << opt::out << std::endl;
    std::cout << "Input rules and regions: " << opt::rules << std::endl;
    std::cout << "Input regions mask file: " << opt::proc_regions << std::endl;
  }

  // set which regions to run
  GRC grv_proc_regions;
  if (opt::proc_regions.length())
    grv_proc_regions.regionFileToGRV(opt::proc_regions);

  // setup the walker
  VariantBamWalker walk(opt::bam, opt::refgenome);
  
  // make the mini rules collection from the rules file
  // this also calls function to parse the BED files
  walk.SetMiniRulesCollection(opt::rules);

  // set the regions to run
  SnowTools::GRC rules_rg = walk.GetMiniRulesCollection().getAllRegions();
  if (grv_proc_regions.size() && rules_rg.size()) // intersect rules regions with mask regions
    rules_rg = rules_rg.intersection(grv_proc_regions);
  else if (grv_proc_regions.size())
    rules_rg = grv_proc_regions; // rules is whole genome, so just make mask instead

  if (rules_rg.size())
    walk.setBamWalkerRegions(rules_rg.asGenomicRegionVector());

  // print out some info
  if (opt::verbose > 0) 
    std::cout << walk << std::endl;;

  // walk the BAM and process
  walk.writeVariantBam();

  // make a bed file
  //if (opt::verbose > 0)
  //  std::cout << "...sending merged regions to BED file" << std::endl;
  //mr->sendToBed("merged_rules.bed");

  //index it
  std::cout << "...writing the index file" << std::endl;
  walk.MakeIndex();

  return;
}

void parseVarOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'v': arg >> opt::verbose; break;
      //case 't': opt::twopass = true; break;
    case 'i': arg >> opt::bam; break;
    case 'o': arg >> opt::out; break;
    case 'g': 
      {
	std::string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;
	if (opt::rules.length())
	  opt::rules += "%";
	opt::rules += "region@" + tmp;
      }
      break;
    case 'c':  // call stats hack
      {
	std::string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;
	if (opt::rules.length())
	  opt::rules += "%";
	opt::rules += "region@" + tmp + ";mate";
      }
      break;

    case 'r': 
      {
	std::string tmp;
	arg >> tmp;
	if (tmp.length() == 0)
	  break;
	if (opt::rules.length())
	  opt::rules += "%";
	else 
	  opt::rules = "region@WG%"; // need to specify a region
	opt::rules += tmp;
      }
      break;
    case 'k': arg >> opt::proc_regions; break;
    case 'f': 
      {
	std::string file;
	arg >> file;
	
	std::string line;
	igzstream iss(file.c_str());
	while(std::getline(iss, line, '\n')) {
	  if (opt::rules.length() && line.length())
	    opt::rules += "%";
	  if (line.length())
	    opt::rules += line;
	}
      }
      break;
    }
  }

  if (opt::bam == "")
    die = true;
  if (opt::out == "")
    die = true;

  // dont stop the run for bad bams for quality checking only
  //opt::perc_limit = opt::qc_only ? 101 : opt::perc_limit;

  // something went wrong, kill
  if (die) {
    std::cout << "\n" << VARIANT_BAM_USAGE_MESSAGE;
    exit(1);
  }

}
