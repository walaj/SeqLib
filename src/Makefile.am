noinst_LIBRARIES = libsnowtools.a

libsnowtools_a_CPPFLAGS = \
        -I.. \
	-I../multifast-v1.4.2 \
        -I../htslib

libsnowtools_a_LIBADD = ../bwa/libbwa.a ../htslib/libhts.a 

libsnowtools_a_SOURCES = GenomicRegion.cpp Histogram.cpp \
			BamWalker.cpp BamStats.cpp MiniRules.cpp gzstream.C SnowUtils.cpp \
			gChain.cpp BamRead.cpp gzstream.C BWAWrapper.cpp BamRead.cpp \
			AlignedContig.cpp STCoverage.cpp DBSnpFilter.cpp PONFilter.cpp \
			BreakPoint.cpp DiscordantCluster.cpp Fractions.cpp \
			../bwa/is.c ../bwa/bwtindex.c \
			../bwa/bwt_gen.c ../bwa/QSufSort.c \
                        ../bwa/bwt.c ../bwa/bwamem_extra.c \
			../multifast-v1.4.2/ahocorasick/ahocorasick.c \
			../multifast-v1.4.2/ahocorasick/node.c 

##bin_PROGRAMS = snowtools

##snowtools_CPPFLAGS = \
##	-I$(top_srcdir)/htslib \
##	-I$(top_srcdir)/multifast-v1.4.2 


## these are the executable versions of snowtools
#snowtools_LDADD = $(top_builddir)/libsnowtools.a 
#snowtools_LDFLAGS = -pthread -std=c++11
#snowtools_SOURCES = SnowTools.cpp CoverageTool.cpp RealignTool.cpp
