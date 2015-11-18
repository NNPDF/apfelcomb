include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDFLAGS	 = 	$(PRJLDFLAGS) ./src/libac_core.a

VPATH=./src
MAIN = appl_comb dis_comb ftdy_comb
DEV = appl_optgrid appl_gridinfo ftdy_hcx cfac_scale

.PHONY: core

all: core $(MAIN)
dev: core $(DEV)

core:
	$(MAKE) -C src/core	

clean:
	-$(RM) -f $(MAIN) $(DEV)
	$(MAKE) clean -C src/core

force_look:
	true