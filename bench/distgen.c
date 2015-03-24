#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>

#define MAXDISTCOUNT 10

int distSize[MAXDISTCOUNT];
int distIter[MAXDISTCOUNT];
int distElems[MAXDISTCOUNT];
double* buffer;
int distsUsed = 0;
int verbose = 0;

void addDist(int size)
{
  int d, dd, mydiff = 0;

  // distances are sorted
  for(d=0; d<distsUsed; d++) {
    if (distSize[d] == size) return; // one dist only once
    if (distSize[d] < size) break;
  }

  assert(distsUsed < MAXDISTCOUNT);
  for(dd = distsUsed; dd >= d; dd--)
    distSize[dd] = distSize[dd-1];

  distSize[d] = size;
  distsUsed++;
}

void initBufs()
{
  int d, i;

  int dElems = (distSize[0]+7)/8;
  buffer = (double*) memalign(64, dElems * 8);
  for(i=0; i<dElems; i++)
    buffer[i] = (double) i;

  for(d=0; d<distsUsed; d++) {
    int dElems = (distSize[d]+7)/8;
    distElems[d] = dElems;
    distIter[d] = distSize[0] / distSize[d];
  }

  if (verbose) {
    fprintf(stderr, "Number of distances: %d\n", distsUsed);
    for(d=0; d<distsUsed; d++)
      fprintf(stderr, "  D%2d: size %d (iter %d)\n",
      	d+1, distSize[d], distIter[d]);
  }
}

void usage(char* argv0)
{
  fprintf(stderr,
	  "Distance Generator\n"
	  "Usage: %s [-<iter>] [<dist1> [<dist2> ... ]]\n\n"
	  "Parameters:\n"
	  "  <iter>       number of times accessing arrays (def. 1000)\n"
	  "  <dist1>, ... different reuse distances (def. 1 dist with 1MB)\n"
	  "Options:\n"
	  "  -h           show this help\n"
	  "  -c           only one access per cache line\n"
	  "  -v           be verbose\n", argv0);
  exit(1);
}

int main(int argc, char* argv[])
{
  int arg, d, i, j, k, dElems;
  int iter = 0;
  int oncePerLine = 0;
  unsigned long aCount = 0;
  double sum = 0.0;

  verbose = 1;  
  for(arg=1; arg<argc; arg++) {
    if (argv[arg][0] == '-') {
      if (argv[arg][1] == 'h') usage(argv[0]);
      if (argv[arg][1] == 'v') { verbose++; continue; }
      if (argv[arg][1] == 'c') { oncePerLine = 1; continue; }
      iter = atoi(argv[arg]+1);
      if (iter == 0) usage(argv[0]);
      continue;
    }
    d = atoi(argv[arg]);
    if (d == 0) usage(argv[0]);
    addDist(d);
  }

  if (distsUsed == 0) addDist(1024*1024);
  if (iter == 0) iter = 1000;

  initBufs();

  if (verbose)
    fprintf(stderr, "Running %d iterations...\n", iter);

  for(i=0; i<iter; i++) {
	sum += buffer[0];
    for(d=0; d<distsUsed; d++) {
    	sum += buffer[0];
		dElems = distElems[d];
		for(k=0; k<distIter[d]; k++) {			
			if (oncePerLine) {
				aCount += dElems/8;
				for(j=0; j<dElems; j+=8)
	  				sum += buffer[j];
			}
	  		else {
	  			aCount += dElems;
	  			for(j=0; j<dElems; j++)
	  				sum += buffer[j];
	  		}
	  	}
    }
  }

  if (verbose)
    fprintf(stderr, "ACount: %lu, sum: %g ...\n", aCount, sum);

  return 0;
}

