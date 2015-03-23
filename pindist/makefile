ifdef PIN_ROOT
CONFIG_ROOT := $(PIN_ROOT)/source/tools/Config
else
CONFIG_ROOT := ../Config
endif
include $(CONFIG_ROOT)/makefile.config
include makefile.rules
include $(TOOLS_ROOT)/Config/makefile.default.rules

#CXXFLAGS=-O3 -fomit-frame-pointer -fno-strict-aliasing -std=c++11
CXXFLAGS=-g -O2 -std=c++11
distcalc: distcalc.cpp
	$(CXX) $(CXXFLAGS) distcalc.cpp -o distcalc
