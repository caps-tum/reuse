/**
 * Distgen
 * Multi-threaded memory access pattern generator
 *
 * Every thread traverses its own nested series of arrays
 * with array sizes as specified. Array sizes correlate
 * to distances in a reuse distance histogram, and fit into
 * given layers (caches) of the memory hierarchy.
 * The idea is to approximate the access pattern of real codes.
 *
 * Copyright 2015 by LRR-TUM
 * Josef Weidendorfer <weidendo@in.tum.de>
 *
 * Licensed under GNU General Public License 2.0 or later.
 * Some rights reserved. See LICENSE
 */

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

struct entry {
  double v;
  struct entry *next;
};

int distSize[MAXDISTCOUNT];
int distIter[MAXDISTCOUNT];
int distBlocks[MAXDISTCOUNT];

// options
int distsUsed = 0;
int verbose = 0;
int iters_perstat = 0;
int tcount = 0;       // number of threads to use
int clockFreq = 2400; // assumed frequency for printing cycles

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
    fprintf(stderr, "  number of distances: %d\n", distsUsed);
    for(d=0; d<distsUsed; d++)
      fprintf(stderr, "    D%2d: size %d (%d traversals per iteration)\n",
      	d+1, distSize[d], distIter[d]);
  }
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

void runBench(struct entry* buffer,
	      int iter, int blocks, int blockDiff,
	      int depChain, int doWrite,
	      double* sum, unsigned long* aCount)
{
  int i, d, k, j, idx, max;
  double lsum, v;
  int idxIncr = blockDiff * BLOCKLEN/sizeof(struct entry);
  int idxMax = blocks * BLOCKLEN/sizeof(struct entry);
  int benchType = depChain + 2*doWrite;

  lsum = *sum;
  for(i=0; i<iter; i++) {
    lsum += buffer[0].v;
    for(d=0; d<distsUsed; d++) {
      lsum += buffer[0].v;
      for(k=0; k<distIter[d]; k++) {
	//fprintf(stderr, "D %d, I %d\n", d, k);
	*aCount += distBlocks[d];
	max = distBlocks[d];
	
	switch(benchType) {
	case 0: // no dep chain, no write
	  idx = 0;
	  for(j=0; j<max; j++) {
	    lsum += buffer[idx].v;
	    idx += idxIncr;
	    if (idx >= idxMax) idx -= idxMax;
	    //fprintf(stderr, " Off %d\n", idx);
	  }
	  break;
	  
	case 1: // dep chain, no write
	  {
	    struct entry* p = buffer;
	    for(j=0; j<max; j++) {
	      lsum += p->v;
	      p = p->next;
	      //fprintf(stderr, " POff %d\n", (int)(p - buffer));
	    }
	  }
	  break;
	  
	case 2: // no dep chain, write
	  idx = 0;
	  for(j=0; j<max; j++) {
	    buffer[idx].v += 1.0;
	    lsum += buffer[idx].v;
	    idx += idxIncr;
	    if (idx >= idxMax) idx -= idxMax;
	    //fprintf(stderr, " Off %d\n", idx);
	  }
	  break;
	  
	case 3: // dep chain, write
	  {
	    struct entry* p = buffer;
	    for(j=0; j<max; j++) {
	      p->v += 1.0;
	      lsum += p->v;
	      p->v = v;
	      p = p->next;
	      //fprintf(stderr, " POff %d\n", (int)(p - buffer));
	    }
	  }
	  break;
	default: assert(0);
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
  if (v > 1000000000000ul)
    sprintf(s,  "%.1f T", 1.0 / 1024.0 / 1024.0 / 1024.0 /1024.0 * v);
  else if (v > 1000000000ul)
    sprintf(s,  "%.1f G", 1.0 / 1024.0 / 1024.0 / 1024.0 * v);
  else if (v > 1000000ul)
    sprintf(s,  "%.1f M", 1.0 / 1024.0 / 1024.0 * v);
  else if (v > 1000ul)
    sprintf(s,  "%.1f K", 1.0 / 1024.0 * v);
  else
    sprintf(s,  "%lu", v);

  return s;
}

int toInt(char* s, int isSize)
{
  char* end;
  double v, f;

  f = isSize ? 1024.0 : 1000.0;
  v = strtof(s, &end);
  if ((*end == 'k') || (*end == 'K'))      v = v * f;
  else if ((*end == 'm') || (*end == 'M')) v = v * f * f;
  else if ((*end == 'g') || (*end == 'G')) v = v * f * f * f;

  return (int) v;
}

void printStats(int ii, double tDiff,
		unsigned long rDiff, unsigned long wDiff)
{
  unsigned long aDiff = rDiff + wDiff;
  double avg = tDiff * tcount / aDiff * 1000000000.0;
  double cTime = 1000.0 / clockFreq;

  fprintf(stderr, " at%5d: ", ii);
  fprintf(stderr,
	  " %5.3fs for %4.1f GB => %5.3f GB/s"
	  " (per core: %6.3f GB/s)\n",
	  tDiff, aDiff * 64.0 / 1000000000.0,
	  aDiff * 64.0 / tDiff / 1000000000.0,
	  aDiff * 64.0 / (tDiff * tcount) / 1000000000.0);
  if (verbose>1)
    fprintf(stderr, "  per access (%lu accesses): %.3f ns (%.1f cycles @ %.1f GHz)\n",
	    aDiff, avg, avg/cTime, 1.0 / 1000.0 * clockFreq);
}

int get_tcount()
{
  static int tc = 0;

  if (tc>0) return tc;

#pragma omp parallel
#ifdef _OPENMP
  tc = omp_get_num_threads();
#else
  tc = 1;
#endif

  return tc;
}

void usage(char* argv0)
{
  fprintf(stderr,
	  "Benchmark with threads accessing their own nested arrays at cache-line granularity\n\n"
	  "Usage: %s [Options] [-<iter>] [<dist1> [<dist2> ... ]]\n"
	  "\nParameters:\n"
	  "  <iter>       number of times (iterations) accessing arrays (def: 1000)\n"
	  "  <dist1>, ... different reuse distances (def: 1 dist with 16MB)\n"
	  "\nOptions:\n"
	  "  -h           show this help\n"
	  "  -p           use pseudo-random access pattern\n"
	  "  -d           traversal by dependency chain\n"
	  "  -w           write after read on traversal\n"
	  "  -c <MHz>     clock frequency in MHz to show cycles per access (def: %d)\n"
	  "  -t <count>   set number of threads to use (def: %d)\n"
	  "  -s <iter>    print perf.stats every few iterations (def: 0 = none)\n"
	  "  -v           be verbose\n", argv0, clockFreq, get_tcount());
  fprintf(stderr,
	  "\nNumbers can end in k/m/g for Kilo/Mega/Giga factor\n");
  exit(1);
}


int main(int argc, char* argv[])
{
  int arg, d, i, j, k, idx;
  int iter = 0, ii;
  int pseudoRandom = 0;
  int depChain = 0;
  int doWrite = 0;
  unsigned long aCount = 0, aCount1;
  int blocks, blockDiff;
  double sum = 0.0;
  int t;
  struct entry* buffer[MAXTHREADS];
  double tt, t1, t2;
  double avg, cTime, gData, gFlops, flopsPA;

  verbose = 0;
  for(arg=1; arg<argc; arg++) {
    if (argv[arg][0] == '-') {
      if (argv[arg][1] == 'h') usage(argv[0]);
      if (argv[arg][1] == 'v') { verbose++; continue; }
      if (argv[arg][1] == 'p') { pseudoRandom = 1; continue; }
      if (argv[arg][1] == 'd') { depChain = 1; continue; }
      if (argv[arg][1] == 'w') { doWrite = 1; continue; }
      if (argv[arg][1] == 'c') {
	if (arg+1<argc) {
	  clockFreq = toInt(argv[arg+1], 0);
	  arg++;
	}
	continue;
      }
      if (argv[arg][1] == 't') {
        if (arg+1<argc) {
          tcount = atoi(argv[arg+1]);
          arg++;
        }
        continue;
      }
      if (argv[arg][1] == 's') {
        if (arg+1<argc) {
	  iters_perstat = toInt(argv[arg+1], 0);
          arg++;
        }
        continue;
      }
      iter = toInt(argv[arg]+1, 0);
      if (iter == 0) {
	fprintf(stderr, "ERROR: expected iteration count, got '%s'\n", argv[arg]+1);
	usage(argv[0]);
      }
      continue;
    }
    d = toInt(argv[arg], 1);
    if (d <= 0) {
      fprintf(stderr, "ERROR: expected distance, got '%s'\n", argv[arg]);
      usage(argv[0]);
    }
    addDist(d);
  }

  if (verbose)
    fprintf(stderr, "Multi-threaded Distance Generator (C) 2015 LRR-TUM\n");
  
  if (distsUsed == 0) addDist(16*1024*1024);
  if (iter == 0) iter = 1000;

  if (tcount == 0) {
    // thread count is the default as given by OpenMP runtime
    tcount = get_tcount();
  }
  else {
    // overwrite thread count of OpenMP runtime
#ifdef _OPENMP
    omp_set_num_threads(tcount);
#else
    // compiled without OpenMP, cannot use more than 1 thread
    if (tcount > 1) {
      fprintf(stderr, "WARNING: OpenMP not available, running sequentially.\n");
      tcount = 1;
    }
#endif
  }
  
  if (iters_perstat == 0) {
    // no intermediate output
    iters_perstat = iter;
  }

  blocks = (distSize[0] + BLOCKLEN - 1) / BLOCKLEN;  
  blockDiff = pseudoRandom ? (blocks * 7/17) : 1;
  blocks = adjustSize(blocks, blockDiff);
  initBufs(blocks);


  // calculate expected number of accesses
  aCount = 0;
  for(d=0; d<distsUsed; d++)
    aCount += distIter[d] * distBlocks[d];

  if (verbose) {
    char sBuf[20], tsBuf[20], acBuf[20], tacBuf[20], tasBuf[20];
    prettyVal(sBuf, BLOCKLEN * blocks);
    prettyVal(tsBuf, BLOCKLEN * blocks * tcount);
    prettyVal(acBuf, aCount);
    prettyVal(tacBuf, aCount * tcount * iter);
    prettyVal(tasBuf, aCount * tcount * iter * 64.0);

    fprintf(stderr, "  buffer size per thread %sB (total %sB), address diff %d\n",
	    sBuf, tsBuf, BLOCKLEN * blockDiff);
    fprintf(stderr, "  accesses per iteration and thread: %s (total %s accs = %sB)\n",
	    acBuf, tacBuf, tasBuf);
  }

  assert(tcount < MAXTHREADS);
  assert(sizeof(struct entry) == 16);
#pragma omp parallel for
  for(t=0; t<tcount; t++) {
    struct entry *next, *buf;
    int idx, blk, nextIdx;
    int idxMax = blocks * BLOCKLEN/sizeof(struct entry);
    int idxIncr = blockDiff * BLOCKLEN/sizeof(struct entry);

    // allocate and initialize used memory
    buffer[t] = (struct entry*) memalign(64, blocks * BLOCKLEN);
    buf = buffer[t];
    for(idx=0; idx < idxMax; idx++) {
      buf[idx].v = (double) idx;
      buf[idx].next = 0;
    }
    // generate dependency chain
    idx = 0;
    for(blk=0; blk < blocks; blk++) {
      nextIdx = idx + idxIncr;
      if (nextIdx >= idxMax) nextIdx -= idxMax;
      //fprintf(stderr, " Blk %d, POff %d\n", blk, nextIdx);
      assert(buf[idx].next == 0);
      buf[idx].next = buf + nextIdx;
      idx = nextIdx;
    }
  }

  fprintf(stderr, "Running %d iterations, %d thread(s) ...\n", iter, tcount);
  if (verbose)
    fprintf(stderr, "  printing statistics every %d iterations\n", iters_perstat);

  aCount = 0;
  tt = wtime();
  t1 = tt;
  ii = 0;
  // loop over chunks of iterations after which statistics are printed
  while(1) {
    aCount1 = aCount;

#pragma omp parallel for reduction(+:sum) reduction(+:aCount)
    for(t=0; t<tcount; t++) {
      double tsum = 0.0;
      unsigned long taCount = 0;

      runBench(buffer[t], iters_perstat, blocks, blockDiff,
	       depChain, doWrite, &tsum, &taCount);

      sum += tsum;
      aCount += taCount;
    }

    t2 = wtime();
    ii += iters_perstat;
    if (ii >= iter) break;

    printStats(ii, t2-t1, aCount - aCount1, doWrite ? (aCount - aCount1) : 0);

    t1 = t2;
  }
  tt = t2 - tt;

  //--------------------------
  // Summary
  //--------------------------

  flopsPA = 1.0;
  if (doWrite) {
    aCount = 2 * aCount;
    flopsPA = .5;
  }

  avg = tt * tcount / aCount * 1000000000.0;
  cTime = 1000.0 / clockFreq;  
  gData = aCount * 64.0 / 1024.0 / 1024.0 / 1024.0;
  gFlops = aCount * flopsPA / 1000000000.0;

  fprintf(stderr, "Summary: throughput %7.3f GB in %.3f s (per core: %.3f GB)\n",
	  gData, tt, gData / tcount);
  fprintf(stderr, "         bandwidth  %7.3f GB/s (per core: %.3f GB/s)\n",
	  gData / tt, gData / tt / tcount);
  fprintf(stderr, "         GFlop/s    %7.3f GF/s (per core: %.3f GF/s)\n",
	  gFlops / tt, gFlops / tt / tcount);
  fprintf(stderr, "         per acc.   %7.3f ns   (%.1f cycles @ %.1f GHz) - \n",
	  avg, avg/cTime, 1.0 / 1000.0 * clockFreq);
  
  if (verbose)
    fprintf(stderr, "         accesses   %lu, sum: %g\n", aCount, sum);


  return 0;
}

