PROGRAM       = acceptanceEff

version       = JTKT
CXX           = g++
CXXFLAGS      = -O -Wall -g -Wno-deprecated -D$(version)
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)
LDFLAGS += -L$(PYTHIA8)/lib -lpythia8
INCS    += -I$(PYTHIA8)/include
CXXFLAGS  += $(INCS)
LDFLAGS += $L -ldl

all:            $(PROGRAM)

$(PROGRAM):    $(PROGRAM).cc
		@echo "Linking $(PROGRAM) ..."
		$(CXX)  -lPhysics -L$(PWD) $(PROGRAM).cc $(CXXFLAGS) $(LDFLAGS) -o $(PROGRAM)
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
