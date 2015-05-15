CXX=g++-4.9

ifdef PIN_ROOT
CONFIG_ROOT := $(PIN_ROOT)/source/tools/Config
else
$(error PIN_ROOT not specified.)
endif
include $(CONFIG_ROOT)/makefile.config
include makefile.rules
include $(TOOLS_ROOT)/Config/makefile.default.rules
#TOOL_CXXFLAGS += -fabi-version=2

installwrapper:
	echo "time \$$PIN_ROOT/pin -t $$PWD/obj-intel64/pinpredict.so \$$*" > pinpredict
	cp pinpredict ~/bin/ && chmod 755 ~/bin/pinpredict


