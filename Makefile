include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) -I./core/
LDFLAGS	 = 	$(PRJLDFLAGS) ./core/libac_core.a

VPATH=./src
MAIN = appl_comb dis_comb ftdy_comb
DEV = appl_optgrid appl_gridinfo ftdy_hcx cfac_scale

.PHONY: core

all: core $(MAIN)
dev: core $(DEV)

core:
	$(MAKE) -C core	

clean:
	-$(RM) -f $(MAIN) $(DEV)
	$(MAKE) clean -C core

force_look:
	true