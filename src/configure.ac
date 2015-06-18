#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(snowtools, 1.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([SnowTools.cpp])
##AC_CONFIG_HEADER([config.h])
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
    boostr=/xchip/gistic/Jeremiah/boost_1_58_0DUMY
    boost_include="-I$boostr"
fi


# Check for the bamtools library pathz
# Bamtools has two different path formats
# If it is built in place the files will be in {path}/lib/ and {path}/include/
# If it is installed (with make install) it will be in {path}/include/bamtools
# and {path}/lib/bamtools. Both variants are handled here
##AC_ARG_WITH(bamtools, AS_HELP_STRING([--with-bamtools=PATH],
##	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]))

##if test "$with_bamtools" -a -d "$with_bamtools"; then
##    bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib -Wl,-rpath,$with_bamtools/lib/bamtools -L$with_bamtools/lib/bamtools"
##    bamtools_include="-I$with_bamtools/include -I$with_bamtools/include/bamtools"
##else
##    bamt=/broad/software/free/Linux/redhat_5_x86_64/pkgs/pezmaster31_bamtools-6708a21
##    bamtools_ldflags="-Wl,-rpath,$bamt/lib -L$bamt/lib -Wl,-rpath,$bamt/lib/bamtools -L$bamt/lib/bamtools"
##    bamtools_include="-I$bamt/include -I$bamt/include/bamtools"
##fi

# Check for the hts library paths
#AC_ARG_WITH(htslib, AS_HELP_STRING([--with-htslib=PATH],
#	[specify directory containing the Heng Li's HTSlib)]))

#if test "$with_htslib" -a -d "$with_htslib"; then
#    htslib_include="-I$with_htslib"
#    htslib_ldflags="-Wl,-rpath,$with_htslib -L$with_htslib"
#else
#    htslibr=	
#    htslib_include="-I$htslibr"
#    htslib_ldflags="-Wl,-rpath,$htslibr -L$htslibr"
#fi

# Check for the aho library paths
#AC_ARG_WITH(aho, AS_HELP_STRING([--with-aho=PATH],
#	[specify directory containing the AhoCorasick string matching algorithm C code)]))

#if test "$with_aho" -a -d "$with_aho"; then
#    aho_include="-I$with_aho"
#    aho_ldflags="-L$with_aho/ahocorasick/build"
#    aho_libs="-lahocorasick"
#else
#    aho=/xchip/gistic/Jeremiah/GIT/isva/ThirdParty/multifast-v1.4.2
#    aho_include="-I$aho"
#    aho_ldflags="-L$aho/ahocorasick"
#fi


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
#AC_SUBST(CPPFLAGS, "$CPPFLAGS $htslib_include $aho_include")
AC_SUBST(CPPFLAGS, "$CPPFLAGS $boost_include")
AC_SUBST(CPPFLAGS, "$CPPFLAGS") 
AC_SUBST(LDFLAGS, "$LDFLAGS")

# We always need to specify to link in bamtools
#AC_SUBST(LIBS, "$LIBS -lhts $aho_libs")
AC_SUBST(LIBS, "$LIBS")

# Make sure the bamtools headers can be found
##AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

# Make sure the ahocorasick headers can be found
#AC_CHECK_HEADERS([ahocorasick/ahocorasick.h],,[AC_MSG_WARN([The Aho-Corasick library must be installed (http://sourceforge.net/projects/multifast/). You can specify its path with the --with-aho=PATH option])])

# Make sure the HTSlib headers can be found
#AC_CHECK_HEADERS([htslib/hts.h],,[AC_MSG_ERROR([The HTSlib library must be installed. You can specify its path with the --with-htslib=PATH option])])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
