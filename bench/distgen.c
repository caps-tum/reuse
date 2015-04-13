#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXDISTCOUNT 10
#define BLOCKLEN 64

int distSize[MAXDISTCOUNT];
int distIter[MAXDISTCOUNT];
int distBlocks[MAXDISTCOUNT];
double* buffer;
int distsUsed = 0;
int verbose = 0;

void addDist(int size)
{
  int d, dd;

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

void initBufs(int blocks)
{
  int d, i;

  for(d=0; d<distsUsed; d++) {
    distBlocks[d] = (distSize[d] + BLOCKLEN - 1) / BLOCKLEN;
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
	  "Usage: %s [Options] [-<iter>] [<dist1> [<dist2> ... ]]\n\n"
	  "Parameters:\n"
	  "  <iter>       number of times accessing arrays (def. 1000)\n"
	  "  <dist1>, ... different reuse distances (def. 1 dist with 1MB)\n"
	  "Options:\n"
	  "  -h           show this help\n"
	  "  -p           use pseudo-random access pattern\n"
	  "  -v           be verbose\n", argv0);
  exit(1);
}

// helper for adjustSize
int gcd(int a, int b)
{
  if (b == 0) return a;
  return gcd(b, a % b);
}

// make sure that gcd(size,diff) is 1 by increasing size, return size
int adjustSize(int size, int diff)
{
  while(gcd(size,diff) >1) size++;
  return size;
}

void runBench(int iter, int blocks, int blockDiff,
	      double* sum, unsigned long* aCount)
{
  int i, d, k, j, idx;
  double* buffer;

  buffer = (double*) memalign(64, blocks * BLOCKLEN);
  for(i=0; i< BLOCKLEN/sizeof(double) * blocks; i++)
    buffer[i] = (double) i;

  for(i=0; i<iter; i++) {
    *sum += buffer[0];
    for(d=0; d<distsUsed; d++) {
      *sum += buffer[0];
      for(k=0; k<distIter[d]; k++) {	
	*aCount += distBlocks[d];
	idx = 0;
	for(j=0; j<distBlocks[d]; j++) {
	  *sum += buffer[idx * BLOCKLEN/sizeof(double)];
	  idx = (idx + blockDiff);
	  if (idx > blocks) idx -= blocks;
	}
      }
    }
  }
}

int main(int argc, char* argv[])
{
  int arg, d, i, j, k, idx;
  int iter = 0;
  int pseudoRandom = 0;
  unsigned long aCount = 0;
  int blocks, blockDiff;
  double sum = 0.0;
  int t, tcount;

  verbose = 1;  
  for(arg=1; arg<argc; arg++) {
    if (argv[arg][0] == '-') {
      if (argv[arg][1] == 'h') usage(argv[0]);
      if (argv[arg][1] == 'v') { verbose++; continue; }
      if (argv[arg][1] == 'p') { pseudoRandom = 1; continue; }
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

  blocks = (distSize[0] + BLOCKLEN - 1) / BLOCKLEN;  
  blockDiff = pseudoRandom ? (blocks * 7/17) : 1;
  blocks = adjustSize(blocks, blockDiff);
  initBufs(blocks);

#pragma omp parallel
#ifdef _OPENMP
  tcount = omp_get_num_threads();
#else
  tcount = 1;
#endif

  if (verbose)
    fprintf(stderr,
	    "Buffer size %d, diff %d. Running %d iterations (%d threads)...\n",
	    BLOCKLEN * blocks, BLOCKLEN * blockDiff, iter, tcount);

#pragma omp parallel for reduction(+:sum) reduction(+:aCount)
  for(t=0; t<tcount; t++) {
    double tsum = 0.0;
    unsigned long taCount = 0;

    runBench(iter, blocks, blockDiff,
	     &tsum, &taCount);

    sum += tsum;
    aCount += taCount;
  }


  if (verbose)
    fprintf(stderr, "ACount: %lu, sum: %g ...\n", aCount, sum);

  return 0;
}

