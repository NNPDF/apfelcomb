################################################################################
#
#	APFELcomb Makefile includes
#
################################################################################

CC= g++
LD= g++

SHELL=/bin/bash

# APPLCOMB paths
RESULTSDIR=   -D  RESULTS_PATH="../nnpdfcpp/data/" -D  DATA_PATH="../nnpdfcpp/data/"
APPLGRIDDIR=  -D  APPL_PATH="../applgrids/"
DATABASEDIR=  -D  DB_PATH="./db/"
ALLDIR= $(RESULTSDIR) $(APPLGRIDDIR) $(DATABASEDIR)

# root
ROOTINCS = $(shell root-config --cflags) 
ROOTLIBS = $(shell root-config --glibs) 

# APFEL
APFELINCS = $(shell apfel-config --cppflags) 
APFELLIBS = $(shell apfel-config --ldflags) 

#LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# applgrid
APPLINCS = -I$(shell applgrid-config --prefix)/include
APPLCLIBS = -L$(shell applgrid-config --prefix)/lib -lAPPLgrid 

# libnnpdf
NNPDFINCLUDE= $(shell nnpdf-config --cppflags)
NNPDFLIBS= $(shell nnpdf-config --ldflags)

# gsl
GSLINCLUDE=$(shell gsl-config --cflags)
GSLLIBS=$(shell gsl-config --libs)

# additional libraries to be included 
PRJLDFLAGS = $(LHAPDFLIBS) $(APPLCLIBS) $(ROOTLIBS) $(APFELLIBS) $(NNPDFLIBS) $(GSLLIBS) -lsqlite3

# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
PRJCXXFLAGS = -Wall -O3 $(ALLDIR) $(LHAPDFINCS) $(APPLINCS) $(ROOTINCS) $(APFELINCS) $(NNPDFINCLUDE) $(GSLINCLUDE)
