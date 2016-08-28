AUTOMAKE_OPTIONS = foreign
SUBDIRS = bwa htslib fermi-lite src

install:  
	mkdir -p bin && mkdir -p include && ar cr bin/libseq.a src/libseqlib.a fermi-lite/libfml.a bwa/libbwa.a htslib/libhts.a

