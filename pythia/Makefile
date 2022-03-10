PROGRAM       = focalPionCorrelation

version       = JTKT
CXX           = /usr/bin/c++
CXXFLAGS      = -O -Wall -g -Wno-deprecated -D$(version)
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)
#CXXFLAGS += $(shell $(FASTJET)/bin/fastjet-config --cxxflags )
#LDFLAGS += $(shell $(FASTJET)/bin/fastjet-config --libs --plugins )
#LDFLAGS += -L$(FASTJET)/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone
LDFLAGS += -L$(PYTHIA8)/lib -lpythia8
INCS    += -I$(PYTHIA8)/include
CXXFLAGS  += $(INCS)
LDFLAGS += $L -ldl

HDRSDICT = include/AliJHMRConst.h include/AliJBaseTrack.h include/AliJHMRHist.h \
include/AliJHMRCorr.h include/AliJHMRPythiaCatalyst.h include/AliJHMRGeantCatalyst.h \
include/AliJHMREvent.h include/AliJHMRCluster.h include/AliJHMRTrack.h

HDRS	+= $(HDRSDICT)  nanoDict.h

SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS) $(PROGRAM).cc
		@echo "Linking $(PROGRAM) ..."
		$(CXX)  -lPhysics -L$(PWD) $(PROGRAM).cc $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)
		@echo "finally done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
	@echo "Compile"
	@echo "$(OUTPUT_OPTION)"
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	@echo "Done"


clean:
		rm -f $(OBJS) core *Dict* $(PROGRAM).o *.d $(PROGRAM) $(PROGRAM).sl

cl:  clean $(PROGRAM)

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcint nanoDict.cc -c -D$(version) $(HDRSDICT) $(INCS)
