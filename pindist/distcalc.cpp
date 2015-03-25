/* Calculate the Stack Reuse Distance Histogram
 * from a trace of block accesses in 'trace.out'
 *
 * (C) 2015, Josef Weidendorfer / LRR-TUM
 * GPLv2+ (see COPYING)
 */

#include <stdio.h>

using namespace std;
#include "dist.cpp"

// must be a power-of-two
#define MEMBLOCKLEN 64

#define CHUNKSIZE 4096
Addr buffer[CHUNKSIZE];

FILE* trace;

void usage(char* argv0)
{
  fprintf(stderr,
	  "Stack Distance Calculation from Trace\n"
	  "Usage: %s [Options] [<trace>]\n\n"
	  "Parameters:\n"
	  "  <trace>      trace file (def. 'trace.out')\n"
	  "Options:\n"
	  "  -h           show this help\n"
	  "  -m <min>     minimal distance\n"
	  "  -s <steps>   number of buckets per doubled distance\n"
	  "  -v           be verbose\n", argv0);
  exit(1);
}

int main(int argc, char* argv[])
{
  int arg;
  int read, ptr;
  unsigned long count;

  int verbose = 0;
  int minDist = 0;
  int doublingSteps = 0;
  const char* tracefile = 0;
  
  for(arg=1; arg<argc; arg++) {
    if (argv[arg][0] == '-') {
      if (argv[arg][1] == 'h') usage(argv[0]);
      if (argv[arg][1] == 'v') { verbose++; continue; }
      if ((argv[arg][1] == 'm') && (arg+1<argc)) {
	minDist = atoi(argv[arg+1]);
	arg++;
	continue;
      }
      if ((argv[arg][1] == 's') && (arg+1<argc)) {
	doublingSteps = atoi(argv[arg+1]);
	arg++;
	continue;
      }

      usage(argv[0]);
    }
    tracefile = argv[arg];
  }

  if (minDist == 0) minDist = 4096;
  if (doublingSteps == 0) doublingSteps = 1;
  if (!trace) tracefile = "trace.out";

  if (verbose)
    fprintf(stderr, "Use trace '%s', min %d, steps %d\n",
	    tracefile, minDist, doublingSteps);

  trace = fopen(tracefile, "r");
  if (trace == 0) {
    fprintf(stderr, "Cannot open '%s'\n", tracefile);
    exit(1);
  }

  double d = minDist;
  double f = pow(2, 1.0/doublingSteps);
  RD_init((int)(d / MEMBLOCKLEN));
  for(d*=f; d< 1024*1024*1024; d*=f)
    RD_addBucket((int)(d / MEMBLOCKLEN));

  count = 0;
  while(1) {
    read = fread(buffer, sizeof(Addr), CHUNKSIZE, trace);
    for(ptr=0; ptr<read; ptr++)
      RD_accessBlock(buffer[ptr]);
    count += read;
    if (read < CHUNKSIZE) break;
  }
  fclose(trace);

  RD_printHistogram(stderr, "", MEMBLOCKLEN);
  fprintf(stderr, "  entries read from trace: %lu\n", count);

  return 0;
}
