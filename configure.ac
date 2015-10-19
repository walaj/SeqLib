s#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(snowtools, 1.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/SnowTools.cpp])
AC_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE([enable])
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
AC_SUBST(CPPFLAGS, "$CPPFLAGS -I. $boost_include")
AC_SUBST(LDFLAGS, "$LDFLAGS")

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS")

# Make sure the bamtools headers can be found
##AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

# Make sure the ahocorasick headers can be found
##AC_CHECK_HEADERS([ahocorasick/ahocorasick.h],,[AC_MSG_WARN([The Aho-Corasick library must be installed (http://sourceforge.net/projects/multifast/). You can specify its path with the --with-aho=PATH option])])

# Make sure the HTSlib headers can be found
##AC_CHECK_HEADERS([htslib/hts.h],,[AC_MSG_ERROR([The HTSlib library must be installed. You can specify its path with the --with-htslib=PATH option])])

# Make sure the HTSlib headers can be found
##AC_CHECK_HEADERS([bwa/bwa.h],,[AC_MSG_ERROR([The BWA library must be installed. You can specify its path with the --with-bwa=PATH option])])

# Make sure the boost headers can be found
##AC_CHECK_HEADERS([boost/icl/interval_set.hpp],,[AC_MSG_WARN([The Boost library must be installed for a few of the interval operations. NOT REQUIRED for VariantBam. Specify its path with the --with-boost=PATH option])])

AC_CONFIG_FILES([Makefile
	         src/Makefile
		 boost/libs/regex/src/Makefile])

AC_OUTPUT
