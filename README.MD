Eer File Format
===============

For the latest eer file format and examples, consult the **EER file format** pdf document in this repository.

How To Build the example TestEerReader
======================================

1. Create a *Build* directory where all the output can be stored. 
2. To generate a 64 bit Visual Studio solution, navigate to the *Build* directory:
    > cmake ../ -DLIBTIFF_INCLUDE_DIR="[Path to Libtiff local include directory]" -DLIBTIFF_LIB_DIR="[Path to Libtiff local lib directory]" -G=[Generator] -A x64

    * For example: 
        > cmake ../ -DLIBTIFF_INCLUDE_DIR="C:\libraries\libtiff\libtiff-4.3.0.1-vc141\include" -DLIBTIFF_LIB_DIR="C:\libraries\libtiff\libtiff-4.3.0.1-vc141\lib\Debug64" -G="Visual Studio 15 2017" -A x64

3. Open the solution in Visual Studio and build the solution. This will create a .lib for the source code (in Build/lib/Debug) and an .exe for the test example (in Build/Debug).

4. To run the TestEerReader.exe use the --help or -h argument on command line to understand what arguments the test needs. 