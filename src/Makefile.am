noinst_LIBRARIES = libsnowtools.a

libsnowtools_a_CPPFLAGS = \
        -I.. \
	-I../multifast-v1.4.2 \
        -I../htslib \
	-I../SGA/Util -I../SGA/Bigraph -I../SGA/SuffixTools \
	-I../SGA/Algorithm -I../SGA/SQG -I../blat/inc -I../blat/jkOwnLib -Wno-write-strings

libsnowtools_a_SOURCES = MiniRules2.cpp jsoncpp.cpp ssw_cpp.cpp ssw.c GenomicRegion.cpp Histogram.cpp RefGenome.cpp SVVCF.cpp \
##			AlignedContig.cpp BreakPoint.cpp \
			BamWalker.cpp BamStats.cpp gzstream.C SnowUtils.cpp \
##			gChain.cpp \
	BamRead.cpp gzstream.C BWAWrapper.cpp BamRead.cpp \
        STCoverage.cpp Fractions.cpp BLATWrapper.cpp \
			../bwa/is.c ../bwa/bwtindex.c \
			../bwa/bwt_gen.c ../bwa/QSufSort.c \
                        ../bwa/bwt.c ../bwa/bwamem_extra.c 

#bin_PROGRAMS = snowtools

#snowtools_CPPFLAGS = \
#	-I$(top_srcdir)/htslib \
#	-I$(top_srcdir)/multifast-v1.4.2 


## these are the executable versions of snowtools
#snowtools_LDADD = libsnowtools.a $(top_srcdir)/htslib/libhts.a $(top_srcdir)/multifast-v1.4.2/ahocorasick/build/libahocorasick.a
##snowtools_LDFLAGS = -std=c++11
#snowtools_SOURCES = snowtools.cpp 
##SnowTools.cpp CoverageTool.cpp RealignTool.cpp
