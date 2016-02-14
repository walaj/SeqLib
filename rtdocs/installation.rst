Installation
------------

The source code is available: https://github.com/jwalabroad/SnowTools

Install the Boost dependency if don't have already:

.. code:: bash
    git clone --recursive https://github.com/boostorg/boost.git
    cd boost
    ./bootstrap.sh --with-libraries=regex,test,filesystem,system
    ./b2

Download and compile the code:

.. code:: bash

    git clone https://github.com/jwalabroad/SnowTools.git
    cd SnowTools
    ./configure --with-boost=<path_to_boost>
    make

C++ Libraries
~~~~~~~~~~~~

To compile HoBBeS, you will need a modern C++ compiler that supports
`c++0x <https://gcc.gnu.org/projects/cxx0x.html>`__ and the dependencies
listed below. I compiled successfully with gcc version 4.9.0 

`GCC, the GNU Compiler <http://gcc.gnu.org>`__

    The GNU Compiler Collection is a compiler system produced by the GNU
    Project supporting various programming languages.

