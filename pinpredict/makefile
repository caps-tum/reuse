ifdef PIN_ROOT
CONFIG_ROOT := $(PIN_ROOT)/source/tools/Config
else
$(error PIN_ROOT not specified.)
endif
include $(CONFIG_ROOT)/makefile.config
include makefile.rules
include $(TOOLS_ROOT)/Config/makefile.default.rules

installwrapper:
	echo "time \$$PIN_ROOT/pin -t $$PWD/obj-intel64/pinpredict.so \$$*" > runpp
	cp runpp ~/bin/ && chmod 755 ~/bin/runpp


