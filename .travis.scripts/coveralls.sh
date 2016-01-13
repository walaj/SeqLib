#!/bin/bash

echo "...trying coveralls"

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
#if [ "$CXX" == "g++" ];
#then
  sudo pip install cpp-coveralls
  cd snow_test
  ./snow_test
  cpp-coveralls -r ../ -e examples -e doxy -e R -e rtdocs -t ${COVERALLS_TOKEN}
#fi
