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

*Provide read/write access to BAM files
*Classes containing individual reads and ways to interact with them (e.g. edit tags)
*Available as an API and a command-line version

SnowTools/BamTools differences
------------------------------
SnowTools/BamTools differences
*Sort/index functionality is independently implemented in BamTools. In SnowTools, the Samtools 
 sort and index functions are called directly.
*BamTools stores quality scores and sequences as strings. In SnowTools, the HTSlib native bam1_t format
 is used instead. This format has a lower memory footprint by using only 4 bits per base, rather than 8. 
 Conversion to C++ style std::string is provided as a function and can be done on the fly.
*BamTools provides the BamMultiReader class for reading multiple BAM files at once, while 
 SnowTools does not currently support this functionality.
*BamTools provides the following CLI tools that are not supported in SnowTools
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
*Jeremiah Wala - Harvard MD-PhD candidate, Bioinformatics and Integrative Genomics
*Steve Schumacher - Research Scientist, Dana Farber Cancer Institute
*Cheng-Zhong Zhang - Research Scientist, Broad Institute
*Marcin Imielinski - Asst Prof, Cornell University

SnowTools vfilter
=================

VariantFilter is a SnowTools module to extract interesting sequencing reads from a BAM file. To save money, 
diskspace and future I/O, one may not want to store an entire BAM on disk after all the relevant VCF, 
MAFs, etc have been created. Instead it would be more efficient to store only those read-pairs 
who intersect some region around the variant locations. Ideally, these regions could be manually inspected
or reanalyzed without having to keep the entire BAM. Having small BAMs is also important in the tool development, where
access to small test files can facilitate a more rapid "build/test" cycle. 

VariantFilter was developed as a one-pass solution for the needs of different NGS tools. For instance, after
initial processing, an indel tool might be interested in storing MAPQ > 0 reads in 30k candidate regions of interest.
A separate SNP tool might find a different 20k regions, while a structural variation tool 
might be interested in only discordant or high-quality clipped reads across the BAM. Alternatively, 
one may want to extract all clipped reads from a BAM while avoiding any reads that lie within LINE/SINE elements. 
VariantFilter uses a series of rules defined on distinct regions in order to handle 
all of these situations with a single pass through the BAM.

In situations where the sequencing or library preparation quality is low, it may be advantageous
to remove poor quality reads before starting the analysis train. VariantFilter handles this by optionally taking into
account Phred base-qualities when making a decision whether to keep a sequencing read. For instance, one might 
only be interested in high quality MAPQ 0 or clipped reads. VariantFilter can be 
setup to apply unique Phred filters to different regions or across the entire genome, all with one-pass.

Note that the capabilities of the [BamTools] command line ``bamtools filter``  
may provide a solution more suitable your needs than VariantFilter. Briefly, ``bamtools filter`` allows you to 
specify a set of filters in JSON format to be applied to a BAM. See the Bamtools documentation for more detail. 

The question then is under what situations would you use ``bamtools filter``, and when would you use VariantFilter? 
Below are a list of scenarios that we feel addresses the different domains of the two tools:

> 1. Extract all MAPQ 0 reads from a BAM - Either tool (prefer ``bamtools filter``)
> 2. Extract all reads in read group A - ``bamtools filter``
> 3. Extract all reads with NM tag >= 4 - Either tool (prefer ``bamtools filter``)
> 4. Extract all reads with NM tag >= 4 in exons - VariantFilter.
> 5. Remove all reads with insert sizes between 100 and 600 bp - VariantFilter
> 6. Extract all reads and mate within 1000bp of a variant or set of genes - VariantFilter
> 7. Extract only high-quality reads with N bases beyong phred score X - VariantFilter
> 8. Reduce a BAM to only high quality reads around your MAFs, VCFs and BED files - VariantFilter

Thus, the main additions of VariantFilter are three-fold:

> 1. The ability to filter specifically on read clipping, orientation and insert size (all important for structural variation), while taking into account the per-base phred quality.
> 2. Use of interval trees to efficiently determine if a read or read mate overlaps a region.
> 3. The ability to provide different rules for different regions, and the ability to provide these regions as common variant files (VCF, MAF, BED)

Example
-------

We ran VariantFilter like this:

```bash
options=(
    --input-bam         big.bam
    --output-bam        small.bam
    --rules-file        rules.vb
    --proc-regions-file small_chr1_mask.bed
)
variant ${options[*]}
```

Syntax
------

This section will describe the syntax used by VariantFilter to specify the cascades of rules and regions 
applied to a BAM. Below is an example of a valid VariantFilter script:

```bash
    ### this is a comment. The line code below defines filters to be applied to each region/rule
    region@WG
    !hardclip;mapped;mapped_mate;isize[0,600];!mapq[10,100]
    !hardclip;mapped;mapped_mate;clip[10,101]
```

##### Region

Let's look first at the ``region`` tag. The region@ keyword marks that what follows is a genomic region, 
which is either the keyword ``WG`` for whole genome, or a VCF, MAF, Callstats or BED file. Regions are 
treated such that they will include any read who overlaps it, even partially. Optionally,
you can specify that your region of interest is a bit bigger than is actually in the file. You can do this by "padding"
the regions around the sites. For example:

``region@myvcf.vcf;pad[1000]``

You can also state that the region applies to reads who don't necessarily overlap the region, but their pair-mate does.

``region@myvcf;pad[1000];mate``

Note that the syntax is such that you must specify the file immediately after the @, following by other options
in any order. 

##### Rules

The next two lines specify a pair of rules, marked with the ``rule@`` tag. 
The default rule is to include every read, and the conditions within the rule are to be  
thought of as exclusion criteria. Note that you can take the complement of a condition 
by prefixing with a ``!``. For example:

```bash
    # do not include hardclipped reads, reads with isize > 600, or reads with mapq between 10 and 100.
    !hardclip;isize[0,600];!mapq[10,100]
    
    # an equivalent specification would be
    !hardclip;mapped;!isize[601,250000000];mapq[0,9]``
```

VariantFilter handles multiple rules in the following way. For each read, VariantFilter 
will cycle through the rules within a region until the read satisfies a rule. When it 
does, it includes the read in the output and stops checking. The logic for the entire collection of 
rules is then as follows:

On a given rule line, the read must satisfy ALL conditions (logical AND)

Across different rules, the read nead only satisfy ONE rule (logical OR)

To illustrate this, note that there is a small discrepancy in the first rule of the above. In the BAM format, 
unmapped reads and reads with unmapped mates are given an insert size of 0. However, in the same rule 
a condition is described to keep all reads with insert sizes 0-600 inclusive. Recalling the AND logic
within a rule, VariantFilter will exclude the read, because it fails the ``mapped`` criteria.

Below is another example which uses the ability of VariantFilter to interpret VCFs and BED files,
and apply rules separately to them.

```bash
    ### declare that region is a VCF file with pads of 1000 on either side of the variant.
    ### use the "mate" keyword to specify that pairs whose mate falls in the region belong to this rule
    region@/home/unix/jwala/myvcf.vcf;mate;pad[1000]
    #### I want to keep all the reads (this the default). Ill be explicit with the "every" keyword
    every
    #### A BED file which gives a list of exons. In here, I just want to keep "variant" reads
    region@/home/unix/jwala/myexonlist.bed 
    ## keep discordant reads
    !isize[0,600];
    ## keep only unmapped reads and their mates
    !mapped;!mapped_mate
    ## or keep if it is hardclipped
    hardclip
    ## keep reads with a mismatch to reference, but with high mapq
    nm[1,101];mapq[30,100]
    
```

##### Global

To reduce redundancy, you can also type a ``global@`` rule anywhere in the stack,
and it will append that rule to everything below. For example, to exclude hardclipped, duplicate, qcfail and 
supplementary reads in every region, you would do:

```bash
    global@!hardclip;!duplicate;!qcfail;!supplementary
    region@WG
    !isize[0,600]
    clip[10,101];mapq[1,60]
    region@myvcf.vcf
```

which is equivalent to

```bash
    region@WG
    !isize[0,600];!hardclip;!duplicate;!qcfail;!supplementary
    clip[10,101];mapq[1,60];!hardclip;!duplicate;!qcfail;!supplementary
    region@myvcf.vcf
    !hardclip;!duplicate;!qcfail;!supplementary
```
	
The global tag will apply through all of the regions. If you want to reset it for everything, just add ``global@every`` 
back onto the stack.

To make things run a little faster, you can set the order so that the more inclusive regions / rules are first. This only
applies if there is an overlap among regions. This is because VariantFilter will move down the list of regions
that apply to this read and stop as soon as it meets an inclusion criteria. I prefer to start with a whole-genome region / rule
set, and then add more fine-mapped regions later.

##### Command Line Script

The usual method of inputing rules is with a VariantFilter script as a text file (passed to
VariantFilter with the ``-r`` flag). However, sometimes it is useful to not have to write an intermediate
file and just feed rules directly in. In that case, just pass a string literal to the -r flag, and VariantFilter
will parse directly. Just separate lines with either a new ``-r`` flag or with a ``%``. For instance, you might run something like the following:

```bash
variant -i big.bam -o small.bam -r 'global@!hardclip' -r 'region@WG%!isize[0,600];%clip[10,101];mapq[1,60]' -r 'region@myvcf.vcf'
```

Note the single quotes so that it is interpreted as a string literal in BASH.

Full list of available rules
----------------------------

```
    #RULE           #EXAMPLE             #DESCRIPTION OF EXAMPLE / FLAG 
    seq	            seq[/home/seqs.txt]  File containing substrings that must be present in the sequence.
    ins             ins[5,101]           Number of inserted bases on the reads (from parsed CIGAR string)
    del             del[10,101]          Number of deleted bases relative to reference (from parsed CIGAR string). 
    nm              nm[0,4]              NM tag from BAM (number of mismatches). e.g. must be 0-4 inclusive
    isize           isize[100,500]       Insert size, where all insert sizes are converted to positive.
    len             len[80,101]          Length of the read following phred trimming
    clip            clip[0,5]            Number of clipped bases following phred trimming
    nbases          nbases[0,5]          Removed reads that have within this range of N bases.
    phred           phred[4,100]         Range of phred scores that are "quality" 
    duplicate       duplicate            Read must be marked as optical duplicate 
    supp            !supp                Read must be primary alignment
    qcfail          !qcfail              Read must note be marked as QC Fail
    fwd_strand      fwd_strand           Read must be mapped to forward strand
    rev_strand      rev_strand           Read must be mapped to reverse strand
    mate_fwd_strand mate_fwd_strand      Mate of read must be mapped to forward strand
    mate_rev_strand mate_rev_strand      Mate of read must be mapped to reverse strand  
    mapped          !mapped              Read must be unmapped
    mapped_mate     mapped_mate          Mate must be mapped
    ff              ff                   Read pair must have forward-forward orientation
    rr              rr                   Read pair must have reverse-reverse orientation
    fr              fr                   Read pair must have forward-reverse orientation (proper)
    rf              rf                   Read pair must have reverse-forward orientation
    ic              ic                   Read pair must have inter-chromosomal mapping
    discordant      discordant[100,600]  Shortcut for !isize[100,600] || rr || ff || rf || ic (!discordant gives "proper" pairs)
```

[license]: https://github.com/broadinstitute/variant-bam/blob/master/LICENSE

[BamTools]: https://raw.githubusercontent.com/wiki/pezmaster31/bamtools/Tutorial_Toolkit_BamTools-1.0.pdf

[API]: http://pezmaster31.github.io/bamtools/annotated.html

