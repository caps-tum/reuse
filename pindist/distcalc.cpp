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

int main()
{
  int read, ptr;
  unsigned long count;

  trace = fopen("trace.out", "r");
  if (trace == 0) {
    fprintf(stderr, "Cannot open trace.out\n");
    exit(1);
  }

  // add buckets [0-1023], [1K - 2K-1], ... [1G - ]
  RD_init(1024 / MEMBLOCKLEN);
  for(int i=2048; i< 1024*1024*1024; i*=2)
    RD_addBucket(i / MEMBLOCKLEN);

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
