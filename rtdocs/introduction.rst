Introduction
------------

HoBBeS is a C++ library designed to give users easy access to the most
powerful and efficient open source methods for sequence assembly, sequence
alignment and BAM/CRAM/SAM acces and processing. HoBBeS is currently
built from a suite of four core methods/libraries:

1. HTSlib_ (https://github.com/samtools/htslib)
2. BWA_ (https://github.com/lh3/bwa)
3. BLAT_ (https://genome.ucsc.edu/cgi-bin/hgBlat?command=start)
4. SGA_ (https://github.com/jts/sga)

Using HoBBeS, a user can read in a series of filtered reads
(e.g. clipped reads only), perform genome assembly, and then do
reference alignment with BWA or BLAT, all with only a few lines of C++
code and with no need to call external programs or write intermediate
files. 

.. _HTSlib https://github.com/samtools/htslib
.. _BWA https://github.com/lh3/bwa
.. _BLAT https://genome.ucsc.edu/cgi-bin/hgBlat?command=start
.. _SGA https://github.com/jts/sga

