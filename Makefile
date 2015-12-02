include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDFLAGS	 = 	$(PRJLDFLAGS) -L ./src -lac_core

MAIN = ./src/appl_comb ./src/dis_comb ./src/ftdy_comb ./src/cfac_scale
DEV = ./src/appl_optgrid ./src/appl_gridinfo ./src/ftdy_hcx 

.PHONY: core

all: core $(MAIN)
	mv ./src/appl_comb ./
	mv ./src/dis_comb ./
	mv ./src/ftdy_comb ./

dev: core $(DEV)
	mv ./src/appl_optgrid ./
	mv ./src/appl_gridinfo ./
	mv ./src/ftdy_hcx ./

core:
	@$(MAKE) -C src/core	

clean:
	-$(RM) -f $(MAIN) $(DEV)
	@$(MAKE) clean -C src/core

force_look:
	true