#CFLAGS=-g -O0
CFLAGS=-O3 -fopenmp -march=native

all: distgen

clean:
	rm -f *.o distgen *~

install: distgen
	cp distgen ~/bin
