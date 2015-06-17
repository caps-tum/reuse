#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXDISTCOUNT 10
#define BLOCKLEN 64
#define MAXTHREADS 64

int distSize[MAXDISTCOUNT];
int distIter[MAXDISTCOUNT];
int distBlocks[MAXDISTCOUNT];
double* buffer;
int distsUsed = 0;
int verbose = 0;

double wtime()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec+1e-6*tv.tv_usec;
}

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
    // each memory block of cacheline size gets accessed
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
	  "Usage: %s [Options] [-<iter>] [<dist1> [<dist2> ... ]]\n"
	  "\nParameters:\n"
	  "  <iter>       number of times accessing arrays (def. 1000)\n"
	  "  <dist1>, ... different reuse distances (def. 1 dist with 1MB)\n"
	  "\nOptions:\n"
	  "  -h           show this help\n"
	  "  -p           use pseudo-random access pattern\n"
	  "  -v           be verbose\n", argv0);
  fprintf(stderr,
	  "\nNumbers can end in k/m/g for Kilo/Mega/Giga factor\n");
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

void runBench(double* buffer, int iter, int blocks, int blockDiff,
	      double* sum, unsigned long* aCount)
{
  int i, d, k, j, idx;

  double lsum = *sum;
  int idxIncr = blockDiff * BLOCKLEN/sizeof(double);
  int idxMax = blocks * BLOCKLEN/sizeof(double);

  for(i=0; i<iter; i++) {
    lsum += buffer[0];
    for(d=0; d<distsUsed; d++) {
      lsum += buffer[0];
      for(k=0; k<distIter[d]; k++) {	
	*aCount += distBlocks[d];
	idx = 0;
	for(j=0; j<distBlocks[d]; j++) {
	  lsum += buffer[idx];
	  idx += idxIncr;
	  if (idx > idxMax) idx -= idxMax;
	}
      }
    }
  }

  *sum = lsum;
}

char* prettyVal(char *s, unsigned long v)
{
  static char str[20];

  if (!s) s = str;
  if (v > 1000000000)
    sprintf(s,  "%.1f G", 1.0 / 1024.0 / 1024.0 / 1024.0 * v);
  else if (v > 1000000)
    sprintf(s,  "%.1f M", 1.0 / 1024.0 / 1024.0 * v);
  else
    sprintf(s,  "%.1f K", 1.0 / 1024.0 * v);

  return s;
}

int toInt(char* s, int isSize)
{
  char* end;
  int d;
  int f = isSize ? 1024 : 1000;

  d = strtol(s, &end, 10);
  if ((*end == 'k') || (*end == 'K')) d = d * f;
  else if ((*end == 'm') || (*end == 'M')) d = d * f * f;
  else if ((*end == 'g') || (*end == 'G')) d = d * f * f * f;
  return d;
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
  double* buffer[MAXTHREADS];
  double tt;

  verbose = 1;  
  for(arg=1; arg<argc; arg++) {
    if (argv[arg][0] == '-') {
      if (argv[arg][1] == 'h') usage(argv[0]);
      if (argv[arg][1] == 'v') { verbose++; continue; }
      if (argv[arg][1] == 'p') { pseudoRandom = 1; continue; }
      iter = toInt(argv[arg]+1, 0);
      if (iter == 0) usage(argv[0]);
      continue;
    }
    d = toInt(argv[arg], 1);
    if (d <= 0) usage(argv[0]);
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

  assert(tcount < MAXTHREADS);
#pragma omp parallel for
  for(t=0; t<tcount; t++) {
    int idx;
    // allocate and initialize used memory
    buffer[t] = (double*) memalign(64, blocks * BLOCKLEN);
    for(idx=0; idx< BLOCKLEN/sizeof(double) * blocks; idx++)
      buffer[t][idx] = (double) idx;
  }

  // calculate expected number of accesses

  aCount = 0;
  for(d=0; d<distsUsed; d++)
    aCount += distIter[d] * distBlocks[d];

  if (verbose) {
    char sBuf[20], tsBuf[20], acBuf[20], tacBuf[20];
    prettyVal(sBuf, BLOCKLEN * blocks);
    prettyVal(tsBuf, BLOCKLEN * blocks * tcount);
    prettyVal(acBuf, aCount);
    prettyVal(tacBuf, aCount * tcount * iter);

    fprintf(stderr, "Buffer size per thread %sB (total %sB), address diff %d\n",
	    sBuf, tsBuf, BLOCKLEN * blockDiff);
    fprintf(stderr, "Accesses per iteration and thread: %s, total %s\n",
	    acBuf, tacBuf);
    fprintf(stderr, "Running %d iterations (%d threads) ...\n",
	    iter, tcount);
  }

  aCount = 0;
  tt = wtime();

#pragma omp parallel for reduction(+:sum) reduction(+:aCount)
  for(t=0; t<tcount; t++) {
    double tsum = 0.0;
    unsigned long taCount = 0;

    runBench(buffer[t], iter, blocks, blockDiff,
	     &tsum, &taCount);

    sum += tsum;
    aCount += taCount;
  }

  tt = wtime() - tt;

  if (verbose) {
    fprintf(stderr, "Finished (ACount: %lu, sum: %g)\n", aCount, sum);
    fprintf(stderr, "Elapsed: %.3fs => %.3f GB/s, %.3f GF/s"
	    " (per core: %.3f GB/s, %.3f GF/s)\n",
	    tt,
	    aCount * 64.0 / tt / 1000000000.0,
	    aCount / tt / 1000000000.0,
	    aCount * 64.0 / (tt * tcount) / 1000000000.0,
	    aCount / (tt * tcount) / 1000000000.0 );
  }

  return 0;
}

