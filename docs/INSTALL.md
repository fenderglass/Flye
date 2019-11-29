Flye Installation
=================

Availability
------------

Flye is available for Linux and MacOS platforms.

Bioconda Releases
-----------------

You can get the latest stable release through Bioconda:

    conda install flye

Alternatively, you can get a release version from the github releases page


Building Requirements
---------------------

* Python 2.7 or 3.5+ (with setuptools package installed)
* C++ compiler with C++11 support (GCC 4.8+ / Clang 3.3+ / Apple Clang 5.0+)
* GNU make
* Git
* Core OS development headers (zlib, ...)


Local building (without installation)
-------------------------------------

You may use the package locally without system installation.
To get and compile the latest git version, run:

    git clone https://github.com/fenderglass/Flye
    cd Flye
    make -j

Then, Flye will be available as:

    python bin/flye


Installing from source
----------------------

To install the Flye package into your system, run:

    git clone https://github.com/fenderglass/Flye
	cd Flye
    python setup.py install

Depending on your OS, you might need to add
```--user``` or ```--prefix``` options to the 
install command for the local installation.

After installation, Flye could be invoked via:

    flye

Optionally, run some tests to ensure that installation was successful:

    python flye/tests/test_toy.py
