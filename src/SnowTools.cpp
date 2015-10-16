#include <string>
#include <getopt.h>
#include <iostream>

//#include "CoverageTool.h"
//#include "VFilterTool.h"
//#include "RealignTool.h"

/*
//static const char* SNOWTOOLS_USAGE_MESSAGE =
static std::string SNOWTOOLS_USAGE_MESSAGE = 
"Usage: snowtools <module> \n\n" +
"Description: Perform a number of operations on BAM files\n" +
    //"  vfilter             Filter a BAM based a series of hierarchical rules\n"
"  coverage            Retrieve binned coverage\n" +
"  realign             Realign reads in single-read mode on BWA-MEM\n" +
"\n";

int main(int argc, char** argv) 
{
  if (argc <= 1) {
    std::cerr << SNOWTOOLS_USAGE_MESSAGE;
    return 0;
  }
  std::string command(argv[1]);
  if (command == "help" || command == "--help" || command == "-h")
    std::cerr << SNOWTOOLS_USAGE_MESSAGE;
  //else if (command == "vfilter")
  //  runVFilter(argc-1, argv+1);
  else if (command == "coverage")
    runCoverage(argc-1, argv+1);
  else if (command == "realign")
    runRealign(argc-1, argv+1);
  else
    std::cerr << std::string(SNOWTOOLS_USAGE_MESSAGE) << "\n\n Input command not recognized" << std::cerr;
  
  return 1;
}

*/
