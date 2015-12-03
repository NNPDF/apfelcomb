include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDFLAGS	 = 	$(PRJLDFLAGS)
LDLIBS   =      ./src/libac_core.a

VPATH=./src
MAIN = appl_comb dis_comb ftdy_comb ./src/cfac_scale
DEV = appl_optgrid appl_gridinfo ftdy_hcx 

.PHONY: core

all: core $(MAIN)
dev: core $(DEV)

core:
	@$(MAKE) -C src/core	

clean:
	-$(RM) -f $(MAIN) $(DEV)
	@$(MAKE) clean -C src/core

force_look:
	true
