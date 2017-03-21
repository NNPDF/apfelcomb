# Core Makefile

include ../../Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) -I ..
LDFLAGS	 = 	$(PRJLDFLAGS)

CORELIB = ../libac_core.a
COREOBJ = fk_appl.o fk_dis.o fk_ftdy.o fk_qcd.o fk_utils.o fk_grids.o

all : $(CORELIB)

$(CORELIB) : $(COREOBJ)
	ar rc $@ $^ && ranlib $@

clean :
	-$(RM) -f $(COREOBJ) $(CORELIB)
force_look :
	true
