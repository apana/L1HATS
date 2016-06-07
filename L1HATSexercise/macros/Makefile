######### GENERAL STUFF: DON NOT CONSIDER
CXX = g++
CXXFLAGS  = -std=c++11
CXXFLAGS  += -Wall -fPIC -O3
#CXXFLAGS  += -g -pg
#CXXFLAGS  += -Wextra
#CXXFLAGS += -Wextra -Weffc++ #Asking for troubles
#CXXFLAGS += -Wno-reorder #Dirty fix
CXXFLAGS  += $(subst -I, -isystem , $(shell root-config --cflags))

LD         = g++
LDFLAGS    =
LIBS       = $(shell root-config --libs)   -lTable

ifneq ($(shell echo $$CMSSW_BASE), )
  CXXFLAGS  += -isystem $(CMSSW_RELEASE_BASE)/src
  CXXFLAGS  += -isystem $(CMSSW_BASE)/src/
  CXXFLAGS  += -isystem $(shell scram tool info boost | awk -F"=" '/INCLUDE=(.*)/{print $$NF}')
  LIBS += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)/ -lFWCoreFWLite -lDataFormatsL1TGlobal
  LIBS += -L$(shell scram tool info boost | awk -F"=" '/LIBDIR=(.*)/{print $$NF}') \
		  -lboost_program_options -lboost_filesystem
endif

PROGRAM = rateExample

OBJS = $(PROGRAM).o $(CMSSW_BASE)/src/L1TriggerDPG/L1Menu/macros/objs/L1Ntuple.o

MKBIN = $(CXX) $(CXXFLAGS) `root-config  --libs --cflags` -lMinuit -lGenVector


$(PROGRAM) : $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	$(LD) $^ $(LIBS)  -o $@
	@echo "done"

all : $(PROGRAM)

$(PROGRAM).cc:

%.o : %.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

%.o : %.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

L1Ntuple.C:


clean:
	rm -f $(PROGRAM) $(OBJS)


