
# Rivet flags
RFLAGS=$(shell rivet-config --cppflags)
RLDFLAGS=$(shell rivet-config --ldflags)
RDATA=$(shell rivet-config --datadir)

# Rivet Library directory - for copying analyses
RLIBS=$(shell rivet-config --libdir)

# APPLgrid flags
AFLAGS=$(shell applgrid-config --cxxflags)
ALDFLAGS=$(shell applgrid-config --ldflags)

# APPLgrid Subprocess directory
SUBDIR=$(shell applgrid-config --share)

# MCgrid flags
MCFLAGS=$(shell pkg-config mcgrid --cflags)
MCLDFLAGS=$(shell pkg-config mcgrid --libs)

# LHAPDF flags
LHAFLAGS= $(shell lhapdf-config --cppflags)
LHALDFLAGS=  $(shell lhapdf-config --ldflags)


#########################################

# Total flags for plugin
PLGFLAGS=$(RFLAGS) $(AFLAGS) $(MCFLAGS)
PLGLDFLAGS=$(RLDFLAGS) $(ALDFLAGS) $(MCLDFLAGS)

# Total flags for test code
TSTFLAGS=$(AFLAGS) $(LHAFLAGS)
TSTLDFLAGS=$(ALDFLAGS) $(LHALDFLAGS)

#########################################

# Plugin
RTARGET = RivetMCgridPlugins.so
RSOURCES = $(wildcard *.cc)

.PHONY: clean install

all:	plugin

plugin: 
	rivet-buildplugin $(RTARGET) $(RSOURCES) $(PLGFLAGS) $(PLGLDFLAGS)
	
clean:
	rm -f *.so 

install:
	cp ./data/* $(RDATA)
	cp ./*.so $(RLIBS)
	cp ./subproc/*.config $(SUBDIR)
