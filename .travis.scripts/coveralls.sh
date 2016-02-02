#!/bin/bash

echo "...running unit tests and code coverage"

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
#if [ "$CXX" == "g++" ];
#then
  sudo pip install cpp-coveralls
  cd snow_test
  ./configure --with-boost=${BOOST_ROOT}
  make
  rm snow_test-snow-test.gcno ## dont get cov on test prog
  ./snow_test
  cpp-coveralls -r ../ -e examples -e doxy -e R -e rtdocs -t ${COVERALLS_TOKEN} -e "snow_test/snow_test.cpp" -e htslib -e bwa -e multifast-v1.4.2
#fi
