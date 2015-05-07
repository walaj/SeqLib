SnowTools
=========

C++ htslib interface and command line tools for interrogating BAM and SAM files.
	
**License:** [GNU GPLv3][license]

Description
-----------

 SnowTools is a C++ package for querying BAM and SAM files, performing 
 BWA-MEM operations in memory, and performing advanced filtering of 
 reads using a hierarchy of rules. At its core, SnowTools is a C++ interface
 to the C htslib project, which provides the code for core BAM/SAM operations. 
 (https://github.com/samtools/htslib)

 SnowTools is built to be extendable. See (...) for examples of how to take advantage of C++
 class extensions to build off of the SnowTools base functionality. 

 SnowTools is available as a stand-alone executable, which wrap a number of tools. 
 However, certain functionality provided in an optimal form by SamTools or 
 BamTools are intentially left out of SnowTools. In cases where improvements or alternative
 approaches are useful, the functionally is partially duplicated (e.g. count).
 
 As noted above, there are many overlaps between this project at the BamTools project from Derek Barnett 
 (https://github.com/pezmaster31/bamtools). As such, SnowTools is provided as a 
 a supplemental library that may be more suited to your individual needs that BamTools, or vice-versa. To
 aid developers in deciding which package is right for them, I have put together a small list of 
 similarities and differences between these two packages

 
SnowTools/BamTools similarity
-----------------------------

> Provide read/write access to BAM files
> Classes containing individual reads and ways to interact with them (e.g. edit tags)
> Available as an API and a command-line version

SnowTools/BamTools differences
------------------------------
SnowTools/BamTools differences
> Sort/index functionality is independently implemented in BamTools. In SnowTools, the Samtools 
 sort and index functions are called directly.
> BamTools stores quality scores and sequences as strings. In SnowTools, the HTSlib native bam1_t format
 is used instead. This format has a lower memory footprint by using only 4 bits per base, rather than 8. 
 Conversion to C++ style std::string is provided as a function and can be done on the fly.
> BamTools provides the BamMultiReader class for reading multiple BAM files at once, while 
 SnowTools does not currently support this functionality.
> BamTools provides the following CLI tools that are not supported in SnowTools
    - convert - Convert between BAM and other formats
    - header - Print the BAM header
    - index - Index a BAM file
    - merge - Merge muliple BAM files
    - random - Select random alignments from BAM file
    - resolve - Resolve paired-end reads IsProperPair flag
    - revert - Remove duplicate marks and restore original base qualities
    - sort - Sorts a BAM file
    - split - Split a BAM file on user-specified property

Installation
------------
```bash
 git clone https://github.com/broadinstitute/isva.git
 cd isva/SnowTools/src
 ./configure
 make
```
 
I have successfully compiled with GCC-4.9.

Support
-------
This code is being developed and maintained by Jeremiah Wala (jwala@broadinstitute.org)

Attributions
------------
> Jeremiah Wala - Harvard MD-PhD candidate, Bioinformatics and Integrative Genomics
> Steve Schumacher - Research Scientist, Dana Farber Cancer Institute
> Cheng-Zhong Zhang - Research Scientist, Broad Institute
> Marcin Imielinski - Asst Prof, Cornell University

