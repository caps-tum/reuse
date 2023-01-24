#CFLAGS=-g -O0
CFLAGS=-O3
UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
	CC=cc
	CFLAGS+=-fopenmp -march=native
endif
ifeq ($(UNAME), Darwin)
	# Apple-provided clang with libomp installed by homebrew
	CC=clang
	CFLAGS+=-Xclang -fopenmp -I/opt/homebrew/opt/libomp/include
	LDFLAGS+=-L/opt/homebrew/opt/libomp/lib -lomp
	#CC=gcc-7
endif

all: distgen

clean:
	rm -f *.o distgen *~

install: distgen
	cp distgen ~/bin
