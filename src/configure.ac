s#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(snowtools, 1.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([SnowTools.cpp])
AC_CONFIG_HEADER([config.h])
##m4_include([m4/m4_ax_openmp.m4])

# Checks for programs.
AC_PROG_CXX ## test for cpp compiler
AC_PROG_CC ## test for C compiler
AC_PROG_RANLIB ## required if libraries are built in package

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Check for libraries
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )

AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
	[specify directory containing the boost library)]))
if test "$with_boost" -a -d "$with_boost"; then
    boost_include="-I$with_boost"
else
    boostr=/xchip/gistic/Jeremiah/boost_1_58_0
    boost_include="-I$boostr"
fi

# Check for the hts library paths
AC_ARG_WITH(htslib, AS_HELP_STRING([--with-htslib=PATH],
	[specify directory containing the Heng Li's HTSlib)]))

if test "$with_htslib" -a -d "$with_htslib"; then
    htslib_include="-I$with_htslib"
    htslib_ldflags="-L $with_htslib"
    ##htslib_ldflags="-Wl,-rpath,$with_htslib -L$with_htslib"
else
    with_htslib="./htslib"
    htslib_include="-I $with_htslib"
    htslib_ldflags="-L $with_htslib"
    ##htslib_ldflags="-Wl,-rpath,$with_htslib -L$with_htslib"
fi
htslib_lib="libhts.a"

# Check for the bwa library paths
#AC_ARG_WITH(bwa, AS_HELP_STRING([--with-bwa=PATH],
#	[dsf]))

#if test "$with_bwa" -a -d "$with_bwa"; then
    #bwa_include="-I$with_bwa"
    #bwa_ldflags="L$with_bwa"
#else
    ##bwa_include="-I $bwar"
#    bwa_ldflags="-L ./bwa"
#fi


# Check for the aho library paths
AC_ARG_WITH(aho, AS_HELP_STRING([--with-aho=PATH],
	[specify directory containing the AhoCorasick string matching algorithm C code)]))

if test "$with_aho" -a -d "$with_aho"; then
    aho_include="-I$with_aho"
    aho_ldflags="-L$with_aho/ahocorasick/build"
    aho_libs="-lahocorasick"
else
    with_aho="./multifast-v1.4.2"
    aho_include="-I $with_aho"
    aho_ldflags="-L $with_aho/ahocorasick/build"
    aho_libs="-lahocorasick"
fi


# Only fail on warnings when the --enable-development flag is passed into configure
AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
##AC_SUBST(AM_CXXFLAGS, "-g -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas -std=c++11")
AC_SUBST(AM_CXXFLAGS, "-g $fail_on_warning -Wno-unknown-pragmas -std=c++11")
AC_SUBST(CXXFLAGS, "-O3")
AC_SUBST(CFLAGS, "-O3")
AC_SUBST(CPPFLAGS, "$CPPFLAGS -I. $htslib_include $aho_include $boost_include")
AC_SUBST(LDFLAGS, "$LDFLAGS -L./bwa $htslib_ldflags $aho_ldflags")

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS $htslib_lib -lahocorasick -lbwa")

# Make sure the bamtools headers can be found
##AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

# Make sure the ahocorasick headers can be found
AC_CHECK_HEADERS([ahocorasick/ahocorasick.h],,[AC_MSG_WARN([The Aho-Corasick library must be installed (http://sourceforge.net/projects/multifast/). You can specify its path with the --with-aho=PATH option])])

# Make sure the HTSlib headers can be found
AC_CHECK_HEADERS([htslib/hts.h],,[AC_MSG_ERROR([The HTSlib library must be installed. You can specify its path with the --with-htslib=PATH option])])

# Make sure the HTSlib headers can be found
AC_CHECK_HEADERS([bwa/bwa.h],,[AC_MSG_ERROR([The BWA library must be installed. You can specify its path with the --with-bwa=PATH option])])

# Make sure the bamtools headers can be found
AC_CHECK_HEADERS([boost/icl/interval_set.hpp],,[AC_MSG_WARN([The Boost library must be installed for a few of the interval operations. NOT REQUIRED for VariantBam. Specify its path with the --with-boost=PATH option])])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
