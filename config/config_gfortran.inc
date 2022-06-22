# Configuration file for compiling software on Ubuntu linux system
# using optimised gfortran compiler

# Set EXEDIR 
EXEDIR = ../exec

# Fortran compiler
F90 ?= gfortran-mp-11
FC = gfortran-mp-11

# Fortran compiler flags
F90FLAGS = -O -cpp
ifeq ($(findstring mpi,$(F90)),mpi)
   F90FLAGS += -DWITH_MPI
endif

FFLAGS = -fdefault-real-8 -fdefault-double-8 -O

FSFLAG = -O

# Debug option
DBG =  

# command to remove files
RM = /bin/rm -f

# extra cleanup command
EXTRA_CLEAN_CMD = :

# Extra files to clean
EXTRA_CLEAN_FILES = *.T

# Module extension
MOD=.mod
