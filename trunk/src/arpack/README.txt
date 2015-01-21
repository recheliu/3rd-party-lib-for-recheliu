This folder is a mirror of ARPACK and ARPACK++ with my patch.

ARPACK:		http://www.caam.rice.edu/software/ARPACK/
ARPACK++:	http://www.ime.unicamp.br/~chico/arpack++/

My patch for ARPACK is to generate 64-bit libraries for visual studio. The patch for ARPACK++ is mainly to fix the compilation errors. Note that ARPACK++ cannot work for complex number since complex in C++ is a template, which is not allowed by C. 

The detail procedure to build ARPACK is in the link below:

http://www.recheliu.org/memo/buildarpackonwindowsforvisualstudio

Here is a recap of the procedure. 

0.	Prerequisites:	CYGWIN 64 and mingw64-x86_64-gcc.

1. 	Modify ARmake.inc in the source code root. I extracted the source code of ARPACK to D:\src\ARPACK. Also I cannot find f77 in the latest MinGW so I used gfortran instead. Thus the following variables should be changed accordingly:
    home = /cygdrive/d/src/ARPACK
    PLAT = x64
    FC = /usr/bin/x86_64-w64-mingw32-gfortran.exe
    FFLAGS    = -O
    RANLIB = /usr/bin/x86_64-w64-mingw32-ranlib.exe

2.	Open a MSYS window and:
    cd /cygdrive/d/src/ARPACK

3.	Change UTIL/second. Remove the sentence: EXTERNAL           ETIME.
    Compile the .f to .o:
    make lib

4.	Wrap the *.o to .dll:
    /usr/bin/x86_64-w64-mingw32-dllwrap.exe --export-all-symbols BLAS/*.o LAPACK/*.o SRC/*.o UTIL/*.o -lm -lgfortran --output-def arpack_x64.def -o arpack_x64.dll
    
5.	Open a Visual Studio (64-bit) Command Prompt and:
    cd d:\src\ARPACK
    
6.	Generate the library:
    lib.exe /machine:X64 /def:arpack_x64.def
	
After building .lib and .dll, adding the directories of related dll. to your path. For CYGWIN 64, the default path should be:

C:\cygwin64\usr\x86_64-w64-mingw32\sys-root\mingw\bin

