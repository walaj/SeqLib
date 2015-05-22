AUTOMAKE_OPTIONS = foreign
SUBDIRS = bwa htslib 

noinst_LIBRARIES = libsnowtools.a

libsnowtools_a_CPPFLAGS = \
	-I$(top_srcdir)/htslib \
	-I$(top_srcdir)/multifast-v1.4.2 

libsnowtools_a_LIBADD = $(top_srcdir)/bwa/libbwa.a $(top_srcdir)/htslib/libhts.a 


libsnowtools_a_SOURCES = GenomicRegion.cpp Histogram.cpp \
			BamWalker.cpp BamStats.cpp MiniRules.cpp gzstream.C SnowUtils.cpp \
			gChain.cpp HTSTools.cpp BamRead.cpp gzstream.C BWAWrapper.cpp \
			bwa/is.c bwa/bwtindex.c \
			bwa/bwt_gen.c bwa/QSufSort.c bwa/bwt.c bwa/bwamem_extra.c \
			multifast-v1.4.2/ahocorasick/ahocorasick.c \
			multifast-v1.4.2/ahocorasick/node.c 

bin_PROGRAMS = snowtools

snowtools_CPPFLAGS = \
	-I$(top_srcdir)/htslib \
	-I$(top_srcdir)/multifast-v1.4.2 

snowtools_LDADD = $(top_builddir)/libsnowtools.a $(top_builddir)/bwa/libbwa.a $(top_builddir)/htslib/libhts.a 
		##$(top_builddir)/multifast-v1.4.2/ahocorasick/build/libahocorasick.a
snowtools_LDFLAGS = -pthread -std=c++11
snowtools_SOURCES = SnowTools.cpp CoverageTool.cpp RealignTool.cpp
