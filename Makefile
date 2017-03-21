include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDLIBS   =  ./src/libac_core.a $(PRJLDFLAGS)

VPATH=./src
MAIN = appl_comb
DEV = grid_analyser
#ftdy_hcx

.PHONY: all dev core clean
	
all: core $(MAIN)
dev: core $(DEV)

core:
	@$(MAKE) -C src/core	

clean:
	-$(RM) -f $(MAIN) $(DEV)
	@$(MAKE) clean -C src/core

force_look:
	true
