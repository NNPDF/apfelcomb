include Makefile.inc

CXXFLAGS = 	$(PRJCXXFLAGS) 
LDLIBS   =  ./src/libac_core.a $(PRJLDFLAGS)

VPATH=./src
MAIN = apfel_comb unnormalization F2_ratio src/cfac_scale 
DEV = appl_optgrid ftdy_hcx

.PHONY: all dev core clean
	
all: core $(MAIN) database
dev: core $(DEV)

core:
	@$(MAKE) -C src/core	

clean:
	-$(RM) -f $(MAIN) $(DEV)
	@$(MAKE) clean -C src/core

database:
	./db/generate_database.sh

force_look:
	true
