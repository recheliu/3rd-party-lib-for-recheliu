###########################################################################
#
#  Program:         ARPACK
#
#  Module:          ARmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 22, 1996
#
#  Modified:
#
#  Send bug reports, comments or suggestions to arpack@caam.rice.edu
#
############################################################################
#
# %---------------------------------%
# |  SECTION 1: PATHS AND LIBRARIES |
# %---------------------------------%
#
#
# %--------------------------------------%
# | You should change the definition of  |
# | home if ARPACK is built some place   | 
# | other than your home directory.      |
# %--------------------------------------%
#
# MOD-BY-LEETEN 2015/01/15:	home = $(HOME)/ARPACK
home = /cygdrive/d/src/ARPACK
# MOD-BY-LEETEN 2015/01/15:	-END
#
#  %--------------------------------------%
#  | The platform identifier to suffix to |
#  | the end of library names             |
#  %--------------------------------------%
#
# MOD-BY-LEETEN 2015/01/15:	# PLAT = SUN4
PLAT = x64
# MOD-BY-LEETEN 2015/01/15-END
	
#
#  %------------------------------------------------------%
#  | The directories to find the various pieces of ARPACK |
#  %------------------------------------------------------%
#
BLASdir      = $(home)/BLAS
LAPACKdir    = $(home)/LAPACK
UTILdir      = $(home)/UTIL
SRCdir       = $(home)/SRC
#
DIRS        = $(BLASdir) $(LAPACKdir) $(UTILdir) $(SRCdir)
#
# %-------------------------------------------------------------------%
# | Comment out the previous line and uncomment the following         |
# | if you already have the BLAS and LAPACK installed on your system. |
# | NOTE: ARPACK assumes the use of LAPACK version 2 codes.           |
# %-------------------------------------------------------------------%
#
#DIRS         = $(UTILdir) $(SRCdir)
#
# %---------------------------------------------------%
# | The name of the libraries to be created/linked to |
# %---------------------------------------------------%
#
ARPACKLIB  = $(home)/libarpack_$(PLAT).a
LAPACKLIB = 
BLASLIB = 
#
ALIBS =  $(ARPACKLIB) $(LAPACKLIB) $(BLASLIB) 
#
# 
# %---------------------------------------------------------%
# |                  SECTION 2: COMPILERS                   |
# |                                                         |
# | The following macros specify compilers, linker/loaders, |
# | the archiver, and their options.  You need to make sure |
# | these are correct for your system.                      |
# %---------------------------------------------------------%
#
#
# %------------------------------%
# | Make our own suffixes' list. |
# %------------------------------%
#
.SUFFIXES:
.SUFFIXES:	.f	.o
#
# %------------------%
# | Default command. |
# %------------------%
#
.DEFAULT:
	@$(ECHO) "Unknown target $@, try:  make help"
#
# %-------------------------------------------%
# |  Command to build .o files from .f files. |
# %-------------------------------------------%
#
.f.o:
	@$(ECHO) Making $@ from $<
	@$(FC) -c $(FFLAGS) $<
#
# %-----------------------------------------%
# | Various compilation programs and flags. |
# | You need to make sure these are correct |
# | for your system.                        |
# %-----------------------------------------%
#
# MOD-BY-LEETEN 2015/01/15-FROM:
# # FC      = f77
# # FFLAGS	= -O -cg89
# MOD-BY-LEETEN 2015/01/15-TO:
FC = /usr/bin/x86_64-w64-mingw32-gfortran.exe
FFLAGS	= -O 
# MOD-BY-LEETEN 2015/01/15-END


LDFLAGS = 
CD      = cd

ECHO    = echo

LN      = ln
LNFLAGS = -s

MAKE    = /bin/make

RM      = rm
RMFLAGS = -f

SHELL   = /bin/sh
#
#  %----------------------------------------------------------------%
#  | The archiver and the flag(s) to use when building an archive   |
#  | (library).  Also the ranlib routine.  If your system has no    |
#  | ranlib, set RANLIB = touch.                                    |
#  %----------------------------------------------------------------%
#
AR = ar 
ARFLAGS = rv
# RANLIB  = touch
# MOD-BY-LEETEN 2015/01/15:	RANLIB   = ranlib
RANLIB = /usr/bin/x86_64-w64-mingw32-ranlib.exe
# MOD-BY-LEETEN 2015/01/15-END

#
# %----------------------------------%
# | This is the general help target. |
# %----------------------------------%
#
help:
	@$(ECHO) "usage: make ?"
