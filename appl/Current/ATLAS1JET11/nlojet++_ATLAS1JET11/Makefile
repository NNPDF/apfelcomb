#Set the Install path as the current directory
# need more standardised directory structure and setup
#INSTPATH=`pwd`
# GPS alternative for when one has a standard nlojet (v4) installation

BASEDIR   = /home/guffanti/software
INSTALLBASE = $(BASEDIR)

NLODIR = $(INSTALLBASE)

INSTPATH= `which nlojet++ | sed s:/bin/nlojet\+\+::`
# CXX = g++
F77 = gfortran
CXXFLAGS += -O3 -fPIC -g 
# CXXFLAGS += -I$(INSTALLBASE)/include/nlo++ -I$(INSTALLBASE)/include
CXXFLAGS += -I$(INSTALLBASE)/include/nlo++ -I$(INSTALLBASE)/include

# HOPPETFLAG = $(shell hoppet.sh $(INSTALLBASE)/include/hoppet_v1.h) -I$(INSTALLBASE)/include

# HOPPETLIB = $(shell hoppet.sh --lib $(INSTALLBASE)/lib/libhoppet_v1.a -L$(INSTALLBASE)/lib -lhoppet_v1) 

HOPPETINCS=$(shell hoppet-config --cxxflags)
HOPPETLIBS=$(shell hoppet-config --libs)
#HOPPETFLAG="-DHOPPET"

CXXFLAGS += $(HOPPETINCS)

FFLAGS += -g 


.PHONY:  clean all
.SUFFIXES : .cxx .o .f .c

.f.o :
	$(F77) $(FFLAGS) -c $<
.cxx.o:	 
	$(CXX) -c $< $(CXXFLAGS) $(ALLMACRO) -I$(INSTALLBASE)/include $(FASTJETINCS) -I$(INSTALLBASE)/include/lhpdf $(APPLCXXFLAGS) 

install : all
all : fillgridP2  gAdd standSimple #standpdf standsubproc # standPowheg standCrossCheck standalphas standalphaspdf


# F2CFLAG = -lg2c
#  -L/usr/lib/gcc/i386-redhat-linux5E/4.1.2/ -lgfortran



# ROOT	
# ROOTINCS = -I$(shell root-config --incdir)

ROOTINCS = -I$(shell root-config --incdir)
ROOTLIBS = $(shell root-config --libs)

CXXFLAGS += $(ROOTINCS)

#LHAPDF
LHAPDFDIR = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

APPLCXXFLAGS = -I$(INSTALLBASE)/include
APPLLIBS     = -L$(INSTALLBASE)/lib -lAPPLgrid

CXXFLAGS += $(APPLCXXFLAGS)

LD_FLAGS += -L$(LHAPDFDIR) -lLHAPDF
LD_FLAGS += $(APPLLIBS)

FASTJET=$(shell fastjet-config --prefix)
FASTJETLIBS = -L$(FASTJET)/lib,-lfastjet
FASTJETINCS = $(shell fastjet-config --cxxflags)

CXXFLAGS += $(FASTJETINCS)

#
#F90 libraries
#
# F90LIB= -L`which lf95 | sed s:bin/lf95:lib:` -lfj9i6 -lfj9f6 -lfj9e6 -lfccx86_6a
#FRTLLIB = $(shell gfortran.sh )
# FRTLLIB = -L/usr/local/gfortran/lib/x86_64  -lgfortran
#FRTLLIB = -lg2c -lgfortran
# FRTLLIB = -lgfortran
F90LIB = $(FRTLLIB)
#PDFCONVLIBS = -L$(INSTPATH)/pdf-conv-1.0/src -ldglap $(F90LIB)
#PDFCONVLIBS = -L$(INSTALLBASE)/lib -ldglap $(F90LIB) $(LHAPDFLIBS)
PDFCONVLIBS = $(HOPPETLIBS) $(F90LIB) $(LHAPDFLIBS)


# macros to setup contact interactions and reference histograms
REF = -DREN_REFERENCE

ALLMACRO  = $(LHA) $(REF)
ALLMACRO1 = $(LHA),$(REF)

#ALL_LIBS = -L$(INSTALLBASE)/lib,-L$(INSTALLBASE)/lib,-lAPPLgrid,-L$(LHAPDFDIR),-lLHAPDF,-L`root-config --libdir`,-lCore,-lCint,-lHist,-lGraf,-lGraf3d,-lGpad,-lTree,-lRint,-lPostscript,-lMatrix,-lPhysics,-lGui,-lRIO,-pthread,-lm,-ldl,-rdynamic
ALL_LIBS = -L$(INSTALLBASE)/lib,-lAPPLgrid,-L$(LHAPDFDIR),-lLHAPDF,`root-config --libs | sed 's/ /,/g'`
#-L`root-config --libdir`,-lCore,-lCint,-lHist,-lGraf,-lGraf3d,-lGpad,-lTree,-lRint,-lPostscript,-lMatrix,-lPhysics,-lGui,-lRIO,-pthread,-lm,-ldl,-rdynamic


#
#   flag for two to tow scattering _P2TO2_
#


NLOLIBTOOLFLAGS = -Wc,$(ARCH),$(ROOTINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/nlo++,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS)

NLOCXXFLAGS3 = -Wx,$(ARCH),$(ALLMACRO1),-D_P2TO3_,-O3,-Wall,$(ROOTINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS),$(FASTJETINCS) 
NLOCXXFLAGS2 = -Wx,$(ARCH),$(ALLMACRO1),-D_P2TO2_,-O3,-Wall,$(ROOTINCS),-I$(INSTALLBASE)/include,-I$(INSTALLBASE)/include/lhpdf,$(APPLCXXFLAGS),$(FASTJETINCS) 
NLOFFLAGS   = -Wf,$(ARCH),-O3,-Wall
NLOLIBFLAGS3 = -Wl,$(ARCH),-lpthread,$(ALL_LIBS),$(FASTJETLIBS),-L$(INSTALLBASE)/lib,-lAPPLgrid,nlogridP3.lo,fjClustering.lo
NLOLIBFLAGS2 = -Wl,$(ARCH),-lpthread,$(ALL_LIBS),$(FASTJETLIBS),-L$(INSTALLBASE)/lib,-lAPPLgrid,nlogridP2.lo,fjClustering.lo


NLOINC = -I$(INSTALLBASE)/include/nlo++ 


.PHONY:  clean all
: all


fjClustering.lo: fjClustering.cpp 
	$(INSTALLBASE)/libexec/nlo++/libtool-nlojet++ --mode=compile $(CXX) $(ALLMACRO) $(CXXFLAGS) -c fjClustering.cpp $(NLOLIBTOOLFLAGS),$(FASTJETINCS)

fillgridP2: LHAPDF.h nlogridP2.lo fjClustering.o fjClustering.lo
	$(INSTALLBASE)/bin/create-nlojet-user lhapdf_string_interface.f fillgrid.cxx -o fillgridP2 $(NLOCXXFLAGS2) $(NLOFFLAGS) $(NLOLIBFLAGS2) 
fillgridP3: LHAPDF.h nlogridP3.lo fjClustering.o fjClustering.lo
	$(INSTALLBASE)/bin/create-nlojet-user lhapdf_string_interface.f fillgrid.cxx -o fillgridP3 $(NLOCXXFLAGS3) $(NLOFFLAGS) $(NLOLIBFLAGS3) 

nlogridP2.lo: nlogrid.cxx 
	$(INSTALLBASE)/libexec/nlo++/libtool-nlojet++ --mode=compile $(CXX) -c nlogrid.cxx -o nlogridP2.lo -D_P2TO2_ $(ALLMACRO) $(NLOLIBTOOLFLAGS)
nlogridP3.lo: nlogrid.cxx 
	$(INSTALLBASE)/libexec/nlo++/libtool-nlojet++ --mode=compile $(CXX) -c nlogrid.cxx -o nlogridP3.lo -D_P2TO3_ $(ALLMACRO) $(NLOLIBTOOLFLAGS)







standRatio: standRatio.o LHAPDF.h
	$(CXX) -o standRatio $(LDFLAGS) $(ALLMACRO)  standRatio.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 
standPowheg: standPowheg.o LHAPDF.h
	$(CXX) -o standPowheg $(LDFLAGS) $(ALLMACRO)  standPowheg.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standScale: standScale.o LHAPDF.h
	$(CXX) -o standScale $(LDFLAGS) $(ALLMACRO)  standScale.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

stand: stand.o LHAPDF.h
	$(CXX) -o $@ $(LDFLAGS) $(ALLMACRO)  $@.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standSimple: standSimple.o LHAPDF.h
	$(CXX) -o $@ $(LDFLAGS) $(ALLMACRO)  $@.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standalphas: standalphas.o LHAPDF.h
	$(CXX) -o $@ $(LDFLAGS) $(ALLMACRO)  $@.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standalphaspdf: standalphaspdf.o LHAPDF.h
	$(CXX) -o $@ $(LDFLAGS) $(ALLMACRO)  $@.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standpdf: standpdf.o LHAPDF.h
	$(CXX) -o $@ $(LDFLAGS) $(ALLMACRO)  $@.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standpdfMRST: standpdfMRST.o LHAPDF.h
	$(CXX) -o standpdfMRST $(LDFLAGS) $(ALLMACRO)  standpdfMRST.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standCrossCheck: standCrossCheck.o LHAPDF.h
	$(CXX) -o standCrossCheck $(LDFLAGS) $(ALLMACRO)  standCrossCheck.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standnnpdf: standnnpdf.o LHAPDF.h
	$(CXX) -o standnnpdf $(LDFLAGS) $(ALLMACRO)  standnnpdf.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standcteq66: standcteq66.o LHAPDF.h
	$(CXX) -o standcteq66 $(LDFLAGS) $(ALLMACRO)  standcteq66.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

standsubproc: standsubproc.o LHAPDF.h
	$(CXX) -o standsubproc $(LDFLAGS) $(ALLMACRO)  standsubproc.o  lhapdf_string_interface.f $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 

gAdd: gAdd.o 
	$(CXX) -o gAdd $(LDFLAGS) $(ALLMACRO)  gAdd.o  $(APPLLIBS) $(ROOTLIBS) $(PDFCONVLIBS) 



clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la standpdf stand standSimple gAdd standsubproc standnnpdf standcteq66  standpdfMRST standRatio *~


