################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: checklibnnpdf rmxeq mkdep mkmain
gen: rmxeq mkdep mkgen
.PHONY: all

GCC = g++ -std=c++11

# APPLCOMB paths
RESULTSDIR=   -D  RESULTS_PATH="../nnpdfcpp/data/" -D  DATA_PATH="../nnpdfcpp/data/"
APPLGRIDDIR=  -D  APPL_PATH="../applgrids/"
DATABASEDIR=  -D  DB_PATH="./"

# main programs and required modules 

MAIN = src/appl_comb src/dis_comb src/ftdy_comb src/appl_optgrid src/appl_gridinfo src/cfac_scale

MODULES = src/fk_appl src/fk_utils src/fk_qcd src/nnpdfdb src/fk_dis src/fk_ftdy

# search path for modules

MDIR = ./

VPATH = .:$(MDIR)/inc:$(MDIR)/

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

GSLINCLUDE=$(shell gsl-config --cflags)
GSLLIBS=$(shell gsl-config --libs)

NNPDFINCLUDE= $(shell nnpdf-config --cppflags)
NNPDFLIBS= $(shell nnpdf-config --ldflags)

# additional include directories
INCPATH = -I./inc $(LHAPDFINCS) $(APPLINCS) $(ROOTINCS) $(APFELINCS) $(NNPDFINCLUDE) $(GSLINCLUDE)

# additional libraries to be included 
LIBS = $(LHAPDFLIBS) $(APPLCLIBS) $(ROOTLIBS) $(RESULTSDIR) $(APPLGRIDDIR) $(APFELLIBS) $(NNPDFLIBS) $(GSLLIBS)

# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS = -Wall -ansi -O3  -lsqlite3 $(RESULTSDIR) $(APPLGRIDDIR) $(DATABASEDIR)

############################## do not change ###################################

SHELL=/bin/bash

CC=$(GCC)

PGMS= $(MAIN) $(DIS) $(MODULES)

INCDIRS = $(INCPATH)

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(LIBS)

-include $(addsuffix .d,$(PGMS))

checklibnnpdf:
	@ if [ $(shell nnpdf-config --safemode) != 1 ]; then  echo "Error: Compile libnnpdf with --enable-safemode"; exit 1; fi;

# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.cc Makefile
	@ $(CC) -MM -ansi $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.cc Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@
	mv $@ ./

$(GEN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@
	mv $@ ./

# produce executables

mkmain: $(MAIN)
mkgen: $(GEN)

# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory 

clean:
	@ -rm -rf ./src/*.d ./src/*.o *~ $(MAIN) $(OPTG) $(DIS) *.in
.PHONY: clean

################################################################################
