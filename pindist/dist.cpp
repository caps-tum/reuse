/* Calculation of a Stack Reuse Distance Histogram
 *
 * (C) 2015, Josef Weidendorfer / LRR-TUM
 * GPLv2+ (see COPYING)
 */

#include "dist.h"

#include <stdio.h>
#include <list>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <algorithm>


// Assertions and consistency check?
#define RD_DEBUG 0

// 2: Huge amount of debug output, 1: checks, 0: silent
#define RD_VERBOSE 0

// use minimal memory block element? Prohibits consistency checks
#define MIN_BLOCKSTRUCT 0

#if MIN_BLOCKSTRUCT
class MemoryBlock {
public:
  MemoryBlock(Addr a) { bucket = 0; } // generated on first access
  void print(char* b)
    { sprintf(b, "block at bucket %d", bucket); }
  void incACount() {}
  unsigned long getACount() { return 1; }

  int bucket;
};
#else
class MemoryBlock {
public:
  MemoryBlock(Addr a)
    { addr = a; bucket = 0; aCount = 1; } // generated on first access
  void print(char* b)
    { sprintf(b, "block %p, bucket %d, aCount %lu", addr, bucket, aCount); }
  void incACount() { aCount++; }
  unsigned long getACount() { return aCount; }

  int bucket;

private:
  Addr addr;
  unsigned long aCount;
};
#endif

list<MemoryBlock> stack;

class Bucket {
public:
  Bucket(int m)  {
    aCount = 0;
    min = m;
    marker = stack.end();
  }
  unsigned long aCount;
  unsigned int min;
  list<MemoryBlock>::iterator marker;
};

vector<Bucket> buckets;
int nextBucket; // when stack is growing, we may enter this bucket

unordered_map<Addr,list<MemoryBlock>::iterator> addrMap;


void RD_init(int min1)
{
  stack.clear();
  addrMap.clear();

  buckets.clear();
  buckets.push_back( Bucket(0) );    // bucket starting with distance 0
  buckets.push_back( Bucket(min1) ); // first real bucket of interest
  buckets.push_back( Bucket(0) );    // for "infinite" distance
  nextBucket = 1;
}

// add distance buckets, starting from smallest (>0)
// only specification of minimal distance required
void RD_addBucket(unsigned int min)
{
  assert(buckets.size() > 2);
  assert(buckets[buckets.size()-2].min < min);
  buckets.insert( buckets.end()-1, Bucket(min) );
}


// do an internal consistency check
void RD_checkConsistency()
{
  unsigned int aCount1, aCount2;
  unsigned int d = 0;
  int b = 0;
  aCount1 = 0;
  aCount2 = buckets[0].aCount;

  if (RD_VERBOSE>0) {
    fprintf(stderr,"\nChecking... (stack size: %lu)\n", stack.size());
    fprintf(stderr,"   START Bucket %d (min depth %u): aCount %lu\n",
	    b, buckets[b].min, buckets[b].aCount );
  }

  list<MemoryBlock>::iterator stackIt;
  for(stackIt = stack.begin(); stackIt != stack.end(); ++stackIt, ++d) {
    if (d == buckets[b+1].min) {
      b++;
      aCount2 += buckets[b].aCount;
      if (RD_VERBOSE>0)
	fprintf(stderr,"   START Bucket %d (min depth %u): aCount %lu\n",
		b, buckets[b].min, buckets[b].aCount );
      assert( stackIt == buckets[b].marker );
    }

    if (RD_VERBOSE>1) {
      static char b[100];
      stackIt->print(b);
      fprintf(stderr,"   %5d: %s\n", d, b);
    }
    aCount1 += stackIt->getACount();
    assert( stackIt->bucket == b );
  }
  assert( nextBucket = b+1 );
  b = buckets.size()-1;
  aCount2 += buckets[b].aCount;
  if (RD_VERBOSE>0) {
    fprintf(stderr,"   Last Bucket %d: aCount %lu\n",
	    b, buckets[b].aCount );
    fprintf(stderr,"   Total aCount: %u\n", aCount1);
  }
#if MIN_BLOCKSTRUCT
  assert( buckets[b].aCount == stack.size() );
#else
  assert( buckets[b].aCount == stack.size() );
  assert( aCount1 == aCount2 );
#endif
}



void moveMarkers(int topBucket)
{
  for(int b=1; b<=topBucket; b++) {
    --buckets[b].marker;
    (buckets[b].marker)->bucket++;
    if (RD_DEBUG)
      assert( (buckets[b].marker)->bucket == b );
  }
}

void handleNewBlock(Addr a)
{
  // new memory block
  stack.push_front( MemoryBlock(a) );
  addrMap[a] = stack.begin();

  if (RD_VERBOSE >1)
    fprintf(stderr," NEW block %p, now %lu blocks seen\n", a, stack.size());

  // move all markers of active buckets
  moveMarkers(nextBucket-1);

  // does another bucket get active?
  if (addrMap.size() <= buckets[nextBucket].min) return;
  if (buckets[nextBucket].min == 0) return; // last bucket reached

  if (RD_VERBOSE >0)
    fprintf(stderr," ACTIVATE bucket %d (next bucket minimum depth %d)\n",
	    nextBucket, buckets[nextBucket+1].min);

  --buckets[nextBucket].marker; // set marker to last entry
  (buckets[nextBucket].marker)->bucket++;
  if (RD_DEBUG)
    assert( (buckets[nextBucket].marker)->bucket == nextBucket );

  nextBucket++;
}

void moveBlockToTop(Addr a, list<MemoryBlock>::iterator blockIt, int bucket)
{
  // move markers of active buckets within range 1 to <bucket>.
  // we need to do this before moving blockIt to top, as
  // there could be a marker on blockIt, which would go wrong
  moveMarkers(bucket);

  // move element pointed to by blockIt to beginning of list
  stack.splice(stack.begin(), stack, blockIt);
  blockIt->incACount();
  blockIt->bucket = 0;

  if (RD_VERBOSE >1)
    fprintf(stderr," MOVED %p from bucket %d to top (aCount %lu)\n",
	    a, bucket, blockIt->getACount());
}


// run stack simulation
// To use a specific block size, ensure that <a> is aligned
void RD_accessBlock(Addr a)
{
  int bucket;
  auto it = addrMap.find(a);
  if (it == addrMap.end()) {
    handleNewBlock(a);
    bucket = buckets.size()-1; // "infinite" distance, put in last bucket
  }
  else {
    // memory block already seen
    auto blockIt = it->second;
    bucket = blockIt->bucket;

    if (blockIt != stack.begin()) {
      moveBlockToTop(a, blockIt, bucket);
    }
    else {
      blockIt->incACount();
      if (RD_VERBOSE >1)
	fprintf(stderr," TOP %p accessed, bucket %d, aCount %lu\n",
		a, bucket, blockIt->getACount());
    }
  }

  buckets[bucket].aCount++;

  if (RD_DEBUG) {
    // run consistency check every 1 million invocations
    static int checkCount = 0;
    checkCount++;
    if (checkCount > 1000000) {
      RD_checkConsistency();
      checkCount = 0;
    }
  }
}


// get statistics
void RD_stat(unsigned long & stack_size, unsigned long & accessCount)
{
  if (RD_DEBUG)
    RD_checkConsistency();

  unsigned long aCount = 0;
  for(const Bucket& b: buckets)
    aCount += b.aCount;

  stack_size = addrMap.size();
  accessCount = aCount;
}


// get resulting histogram
// Repeatly call RD_get_hist, start with bucket 0.
// Returns next bucket or 0 if this was last
int RD_get_hist(unsigned int b,
		unsigned int & min, unsigned long & accessCount)
{
  if (RD_DEBUG)
    RD_checkConsistency();

  assert((b>=0) && (b < buckets.size()));
  min = buckets[b].min;
  accessCount = buckets[b].aCount;
  if (b == buckets.size()-1) return 0;
  return b+1;
}
