noinst_LIBRARIES = libsnowtools.a

libsnowtools_a_CPPFLAGS = \
        -I$(top_srcdir) \
	-I$(top_srcdir)/multifast-v1.4.2 \
        -I$(top_srcdir)/htslib

libsnowtools_a_LIBADD = $(top_srcdir)/bwa/libbwa.a $(top_srcdir)/htslib/libhts.a 

libsnowtools_a_SOURCES = GenomicRegion.cpp Histogram.cpp \
			BamWalker.cpp BamStats.cpp MiniRules.cpp gzstream.C SnowUtils.cpp \
			gChain.cpp BamRead.cpp gzstream.C BWAWrapper.cpp BamRead.cpp \
			AlignedContig.cpp STCoverage.cpp DBSnpFilter.cpp \
			BreakPoint.cpp DiscordantCluster.cpp \
			$(top_srcdir)/bwa/is.c $(top_srcdir)/bwa/bwtindex.c \
			$(top_srcdir)/bwa/bwt_gen.c $(top_srcdir)/bwa/QSufSort.c \
                        $(top_srcdir)/bwa/bwt.c $(top_srcdir)/bwa/bwamem_extra.c \
			$(top_srcdir)/multifast-v1.4.2/ahocorasick/ahocorasick.c \
			$(top_srcdir)/multifast-v1.4.2/ahocorasick/node.c 

##bin_PROGRAMS = snowtools

##snowtools_CPPFLAGS = \
##	-I$(top_srcdir)/htslib \
##	-I$(top_srcdir)/multifast-v1.4.2 


## these are the executable versions of snowtools
#snowtools_LDADD = $(top_builddir)/libsnowtools.a 
#snowtools_LDFLAGS = -pthread -std=c++11
#snowtools_SOURCES = SnowTools.cpp CoverageTool.cpp RealignTool.cpp
