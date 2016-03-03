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
 * Copyright 2015-2016 by LRR-TUM
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


#define VERSION "0.1"
#define MAXDISTCOUNT 32
#define BLOCKLEN 64
#define MAXTHREADS 256

typedef unsigned long long u64;

struct entry {
  double v;
  struct entry *next;
};

u64 distSize[MAXDISTCOUNT];
u64 distBlocks[MAXDISTCOUNT];
int distIter[MAXDISTCOUNT];

// options (to be reset to default if 0)
int distsUsed = 0;
int verbose = 0;
int iters_perstat = 0;
int tcount = 0;       // number of threads to use
u64 clockFreq = 0;    // assumed frequency for printing cycles
int pseudoRandom = 0;
int depChain = 0;
int doWrite = 0;
int iter = 0;

// defaults
char* clockFreqDef = "2.4G";

double wtime()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec+1e-6*tv.tv_usec;
}

void addDist(u64 size)
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

void initBufs()
{
  int d;

  for(d=0; d<distsUsed; d++) {
    // each memory block of cacheline size gets accessed
    distBlocks[d] = (distSize[d] + BLOCKLEN - 1) / BLOCKLEN;
    distIter[d] = (int) (distSize[0] / distSize[d]);
  }

  if (verbose) {
    fprintf(stderr, "  number of distances: %d\n", distsUsed);
    for(d=0; d<distsUsed; d++)
      fprintf(stderr, "    D%2d: size %llu (%d traversals per iteration)\n",
      	d+1, distSize[d], distIter[d]);
  }
}


// helper for adjustSize
int gcd(u64 a, u64 b)
{
  if (b == 0) return a;
  return gcd(b, a % b);
}

// make sure that gcd(size,diff) is 1 by increasing size, return size
int adjustSize(u64 size, u64 diff)
{
  while(gcd(size,diff) >1) size++;
  return size;
}

void runBench(struct entry* buffer,
	      int iter, u64 blocks, u64 blockDiff,
	      int depChain, int doWrite,
	      double* sum, u64* aCount)
{
  int i, d, k;
  u64 j, idx, max;
  double lsum, v;
  u64 idxIncr = blockDiff * BLOCKLEN/sizeof(struct entry);
  u64 idxMax = blocks * BLOCKLEN/sizeof(struct entry);
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

char* prettyVal(char *s, u64 v)
{
  static char str[50];

  if (!s) s = str;
  if (v >= 1000000000000ull)
    sprintf(s,  "%.1f T", 1.0 / 1000000.0 / 1000000.0 * v);
  else if (v >= 1000000000ull)
    sprintf(s,  "%.1f G", 1.0 / 1000000000.0 * v);
  else if (v >= 1000000ull)
    sprintf(s,  "%.1f M", 1.0 / 1000000.0 * v);
  else if (v >= 1000ull)
    sprintf(s,  "%.1f K", 1.0 / 1000.0 * v);
  else
    sprintf(s,  "%llu", v);

  return s;
}

u64 toU64(char* s, char def)
{
  u64 num = 0, denom = 1;
  char fac = def;

  while((*s >= '0') && (*s <='9')) {
    num = 10*num + (*s - '0');
    s++;
  }
  if (*s == '.') {
    s++;
    while((*s >= '0') && (*s <='9')) {
      num = 10*num + (*s - '0');
      denom = 10*denom;
      s++;
    }
  }

  if (*s) fac = *s;
  if ((fac == 'k') || (fac == 'K'))      num = num * 1000;
  else if ((fac == 'm') || (fac == 'M')) num = num * 1000000;
  else if ((fac == 'g') || (fac == 'G')) num = num * 1000000000ull;
  num = num / denom;

  return num;
}

void printStats(int ii, double tDiff, u64 rDiff, u64 wDiff)
{
  u64 aDiff = rDiff + wDiff;
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
    fprintf(stderr, "  per access (%llu accesses): %.3f ns (%.1f cycles @ %.1f GHz)\n",
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
	  "  <dist1>, ... reuse distances in bytes (without unit assume MiB)\n"
	  "               (default: 1 distance with 16 MiB)\n"
	  "\nOptions:\n"
	  "  -h           show this help\n"
	  "  -r           use (pseudo-)random access pattern\n"
	  "  -d           traversal by dependency chain\n"
	  "  -w           write after read on each access\n"
	  "  -c <freq>    clock frequency in Hz to show cycles per access (def: %s)\n"
	  "  -t <count>   set number of threads to use (def: %d)\n"
	  "  -s <iter>    print perf.stats every few iterations (def: 0 = none)\n"
	  "  -v           be verbose\n", argv0, clockFreqDef, get_tcount());
  fprintf(stderr,
	  "\nNumbers can end in k/m/g for Kilo/Mega/Giga factor\n");
  exit(1);
}

// this sets global options
void parseOptions(char argc, char* argv[])
{
  int arg, pos;
  u64 dist;

  for(arg=1; arg<argc; arg++) {
    if (argv[arg][0] == '-') {
      for(pos=1; argv[arg][pos]; pos++) {
        if (argv[arg][pos] == 'h') usage(argv[0]);
        if (argv[arg][pos] == 'v') { verbose++; continue; }
        if (argv[arg][pos] == 'r') { pseudoRandom = 1; continue; }
        if (argv[arg][pos] == 'd') { depChain = 1; continue; }
        if (argv[arg][pos] == 'w') { doWrite = 1; continue; }

        if (argv[arg][pos] == 'c') {
	  if (argv[arg][pos+1])
	    clockFreq = toU64(argv[arg]+pos+1,'g');
	  else if (arg+1<argc) {
	    clockFreq = toU64(argv[arg+1],'g');
	    arg++;
	  }
        }
        else if (argv[arg][pos] == 't') {
	  if (argv[arg][pos+1])
	    tcount = atoi(argv[arg]+pos+1);
	  else if (arg+1<argc) {
            tcount = atoi(argv[arg+1]);
            arg++;
          }
        }
        else if (argv[arg][pos] == 's') {
	  if (argv[arg][pos+1])
	    iters_perstat = (int) toU64(argv[arg]+pos+1, 0);
	  else if (arg+1<argc) {
	    iters_perstat = (int) toU64(argv[arg+1], 0);
            arg++;
          }
        }
	else {
          iter = (int) toU64(argv[arg]+pos, 0);
          if (iter == 0) {
	    fprintf(stderr, "ERROR: expected iteration count, got '%s'\n",
		    argv[arg]+pos);
	    usage(argv[0]);
          }
	}
	break; // next arg
      }
      continue;
    }

    dist = toU64(argv[arg], 'm');
    if (dist == 0) {
      fprintf(stderr, "ERROR: expected distance, got '%s'\n", argv[arg]);
      usage(argv[0]);
    }
    addDist(dist);
  }

  // set to defaults if values were not provided
  
  if (distsUsed == 0) addDist(16000000);
  if (iter == 0) iter = 1000;
  if (clockFreq == 0)
    clockFreq = toU64(clockFreqDef,'g');

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


}

int main(int argc, char* argv[])
{
  int i, j, k, d;
  int ii;
  u64 aCount = 0, aCount1;
  u64 blocks, blockDiff;
  double sum = 0.0;
  int t;
  struct entry* buffer[MAXTHREADS];
  double tt, t1, t2;
  double avg, cTime, gData, gFlops, flopsPA;

  parseOptions(argc, argv);

  if (verbose) {
    fprintf(stderr, "Multi-threaded Distance Generator ");
    fprintf(stderr, "(C) 2015-2016 LRR-TUM v%s\n", VERSION);
  }
  blocks = (distSize[0] + BLOCKLEN - 1) / BLOCKLEN;  
  blockDiff = pseudoRandom ? (blocks * 7/17) : 1;
  blocks = adjustSize(blocks, blockDiff);
  initBufs();

  // calculate expected number of accesses
  aCount = 0;
  for(d=0; d<distsUsed; d++)
    aCount += distIter[d] * distBlocks[d];
  if (doWrite)
    aCount += aCount;

  if (verbose) {
    char sBuf[20], tsBuf[20], acBuf[20], tacBuf[20], tasBuf[20];
    prettyVal(sBuf, BLOCKLEN * blocks);
    prettyVal(tsBuf, BLOCKLEN * blocks * tcount);
    prettyVal(acBuf, aCount);
    prettyVal(tacBuf, aCount * tcount * iter);
    prettyVal(tasBuf, aCount * tcount * iter * 64.0);

    fprintf(stderr, "  buffer size per thread %siB (total %siB), address diff %llu\n",
	    sBuf, tsBuf, BLOCKLEN * blockDiff);
    fprintf(stderr, "  accesses per iteration and thread: %s (total %s accs = %siB)\n",
	    acBuf, tacBuf, tasBuf);
  }

  assert(tcount < MAXTHREADS);
  assert(sizeof(struct entry) == 16);

  //--------------------------
  // Initialization
  //--------------------------

#pragma omp parallel for
  for(t=0; t<tcount; t++) {
    struct entry *next, *buf;
    u64 idx, blk, nextIdx;
    u64 idxMax = blocks * BLOCKLEN/sizeof(struct entry);
    u64 idxIncr = blockDiff * BLOCKLEN/sizeof(struct entry);

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

  //--------------------------
  // Run benchmark
  //--------------------------

  fprintf(stderr, "Do %d %s%s%siters, %d thr(s) [",
	  iter,
	  pseudoRandom ? "random " : "",
	  depChain ? "chained " : "",
	  doWrite ? "writing " : "",
	  tcount);
  for(d=0; d<distsUsed; d++)
    fprintf(stderr, "%s%s", (d==0) ? "":", ", prettyVal(0, distSize[d]));
  fprintf(stderr, "] (total %siB) ...\n",
		  prettyVal(0, BLOCKLEN * blocks * tcount));

  if (verbose)
    fprintf(stderr, "  printing statistics every %d iterations\n",
	    iters_perstat);

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
      u64 taCount = 0;

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
  cTime = 1000000000.0 / clockFreq;
  gData = aCount * 64.0 / 1000000000.0;
  gFlops = aCount * flopsPA / 1000000000.0;

  fprintf(stderr, "Summary: throughput %7.3f GiB in %.3f s (per core: %.3f GiB)\n",
	  gData, tt, gData / tcount);
  fprintf(stderr, "         bandwidth  %7.3f GiB/s (per core: %.3f GiB/s)\n",
	  gData / tt, gData / tt / tcount);
  fprintf(stderr, "         GFlop/s    %7.3f GF/s (per core: %.3f GF/s)\n",
	  gFlops / tt, gFlops / tt / tcount);
  fprintf(stderr, "         per acc.   %7.3f ns   (%.1f cycles @ %.2f GHz)\n",
	  avg, avg/cTime, 1.0 / 1000000000.0 * clockFreq);
  
  if (verbose)
    fprintf(stderr, "         accesses   %llu, sum: %g\n", aCount, sum);

  return 0;
}

