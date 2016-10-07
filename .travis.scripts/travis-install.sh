#!/usr/bin/env bash

set -e
set -o pipefail

## print info aobut compiler
TC=`$COMPILER --version`
echo "True compiler is $TC"

## only build for one compiler
#VALID=`g++ --version | grep 4.9 | wc -l`
#if [[ "$VALID" -eq "0" ]]; then
#  exit 0;
#fi

if [ -d "${BOOST_ROOT}" ]; then
  (cd "${BOOST_ROOT}"
    ./bootstrap.sh --with-libraries="${BOOST_LIBS}"
    ./b2 threading=multi --prefix="${BOOST_ROOT}" -d0 install
  )
fi
