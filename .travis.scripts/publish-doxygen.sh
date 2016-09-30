#!/bin/bash

## only build for one compiler
if [ "$COMPILER" == "g++-4.9" ] && [ "$TRAVIS_BRANCH" == "master" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ]; 
then 
    
    echo -e "Downloading latest Doxygen...";
    cd ${HOME};
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/travis/doxygen_1.8.8-1_amd64.deb
    sudo dpkg --install doxygen_1.8.8-1_amd64.deb
    cd ${HOME}/build/walaj/SeqLib;
    doxygen
    
    echo -e "Publishing doxygen...\n";
    git config --global user.email "travis@travis-ci.org";
    git config --global user.name "travis-ci";
    git clone --branch=gh-pages https://${GH_TOKEN}@github.com/walaj/SeqLib gh-pages;
    cd gh-pages;
    rm -rf doxygen/;
    mv ../docs/html doxygen/;
    echo "...listing doc" 
    p=`pwd`
    echo "...current location ${pwd}"
    cd doxygen && ls && cd ..
    git add doxygen/;
    git commit -am "Latest doxygen documentation on successful travis build ${TRAVIS_BUILD_NUMBER} auto-pushed";
    git push origin gh-pages 
    
    echo -e "Published doxygen.\n"
    
fi
#fi
