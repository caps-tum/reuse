#include <stdio.h>
#include "pin.H"
#include <list>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <algorithm>

// Assertions and consistency check?
#define DEBUG 1

// 2: Huge amount of debug output, 1: checks, 0: silent
#define VERBOSE 1

// uses INS_IsStackRead/Write: misleading with -fomit-frame-pointer
#define IGNORE_STACK 1

// collect addresses in chunk buffer before?
#define MERGE_CHUNK 0
#define CHUNKSIZE 4096

// must be a power-of-two
#define MEMBLOCKLEN 64

typedef void* Addr;

typedef struct
{
  Addr addr;
  int bucket;
  unsigned int aCount;
} MemoryBlock;

list<MemoryBlock> stack;

class Bucket {
public:
  Bucket(int minDist)  {
    aCount = 0;
    minDepth = minDist / MEMBLOCKLEN;
    marker = stack.end();
  }
  int minDist() const { return minDepth * MEMBLOCKLEN; }

  unsigned long aCount;
  unsigned int minDepth;
  list<MemoryBlock>::iterator marker;
};

vector<Bucket> buckets;
int nextBucket; // when stack is growing, we may enter this bucket

unordered_map<Addr,list<MemoryBlock>::iterator> addrMap;

// used with IGNORE_STACK
unsigned long stackAccesses;

/* ===================================================================== */
/* Handle Memory block access (aligned at multiple of MEMBLOCKLEN)       */
/* ===================================================================== */

void accessBlock(Addr a)
{
  int bucket;
  auto it = addrMap.find(a);
  if (it == addrMap.end()) {
    // new memory block
    stack.push_front( {a,0,1} );
    addrMap[a] = stack.begin();
    bucket = buckets.size()-1; // "infinite" distance, put in last bucket

    if (VERBOSE >1)
      fprintf(stderr," NEW N block %p, Bucket %d\n", a, bucket);

    // move all markers of active buckets
    for(int b=1; b<nextBucket; b++) {
      --buckets[b].marker;
      (buckets[b].marker)->bucket++;
      if (DEBUG)
	assert( (buckets[b].marker)->bucket == b );
    }

    // does another bucket get active?
    if (stack.size() > buckets[nextBucket].minDepth) {
      if (VERBOSE >0)
	fprintf(stderr," MARK bucket %d (next bucket minDepth %d)\n",
		nextBucket, buckets[nextBucket+1].minDepth);
      --buckets[nextBucket].marker; // set marker to last entry
      (buckets[nextBucket].marker)->bucket++;
      if (DEBUG)
	assert( (buckets[nextBucket].marker)->bucket == nextBucket );
      nextBucket++;
    }
  }
  else {
    // memory block already seen
    auto blockIt = it->second;
    bucket = blockIt->bucket;
    if (blockIt != stack.begin()) {
      // move markers of active buckets within range 1 to <bucket>.
      // we need to do this before moving blockIt to top, as
      // there could be a marker on blockIt, which would get invalid
      for(int b=1; b<=bucket; b++) {
	--buckets[b].marker;
	(buckets[b].marker)->bucket++;
	if (DEBUG)
	  assert( (buckets[b].marker)->bucket == b );
      }

      unsigned int aCount = blockIt->aCount +1;
      stack.erase(blockIt);
      stack.push_front( {a,0,aCount} );
      addrMap[a] = stack.begin();
      if (VERBOSE >1)
	fprintf(stderr," MOVED %p, Bucket %d, aCount %u\n",
		a, bucket, aCount);
    }
    else {
      blockIt->aCount++;
      if (VERBOSE >1)
	fprintf(stderr," TOP %p accessed, Bucket %d, aCount %u\n",
		a, bucket, blockIt->aCount);
    }
  }

  buckets[bucket].aCount++;
}

void checkStack()
{
  unsigned int aCount1, aCount2;
  unsigned int d = 0;
  int b = 0;
  aCount1 = 0;
  aCount2 = buckets[0].aCount;

  if (VERBOSE>0) {
    fprintf(stderr,"\nChecking... (stack size: %lu)\n", stack.size());
    fprintf(stderr,"   START Bucket %d (minDepth %u): aCount %lu\n",
	    b, buckets[b].minDepth, buckets[b].aCount );
  }

  list<MemoryBlock>::iterator stackIt;
  for(stackIt = stack.begin(); stackIt != stack.end(); ++stackIt, ++d) {
    if (d == buckets[b+1].minDepth) {
      b++;
      aCount2 += buckets[b].aCount;
      if (VERBOSE>0)
	fprintf(stderr,"   START Bucket %d (minDepth %u): aCount %lu\n",
		b, buckets[b].minDepth, buckets[b].aCount );
      assert( stackIt == buckets[b].marker );      
    }

    if (VERBOSE>1)
      fprintf(stderr,"   %5d: Block %p, Bucket %d, aCount %u\n",
	      d, stackIt->addr, stackIt->bucket, stackIt->aCount);
    aCount1 += stackIt->aCount;
    assert( stackIt->bucket == b );
  }
  assert( nextBucket = b+1 );
  b = buckets.size()-1;
  aCount2 += buckets[b].aCount;
  if (VERBOSE>0) {
    fprintf(stderr,"   Last Bucket %d: aCount %lu\n",
	    b, buckets[b].aCount );
    fprintf(stderr,"   Total aCount: %u\n", aCount1);
    fprintf(stderr,"   Ignored stack accesses: %lu\n\n", stackAccesses);
  }
  assert( buckets[b].aCount == stack.size() );
  assert( aCount1 == aCount2 );
}

#if MERGE_CHUNK
void accessMerging(Addr a)
{
  static Addr mergeBuffer[CHUNKSIZE];
  static int ptr = 0;

  if (ptr < CHUNKSIZE) {
    mergeBuffer[ptr++] = a;
    return;
  }
  sort(mergeBuffer,mergeBuffer+CHUNKSIZE);
  for(ptr=0; ptr<CHUNKSIZE; ptr++) {
    accessBlock(mergeBuffer[ptr]);
  }
  ptr = 0;
}
#define accessBlock accessMerging
#endif


/* ===================================================================== */
/* Direct Callbacks                                                      */
/* ===================================================================== */

void memAccess(ADDRINT addr, UINT32 size)
{
  static int checkCount = 0;

  Addr a1 = (void*) (addr & ~(MEMBLOCKLEN-1));
  Addr a2 = (void*) ((addr+size-1) & ~(MEMBLOCKLEN-1));
  if (a1 == a2) {
    if (VERBOSE >1)
      fprintf(stderr," => %p\n", a1);
    accessBlock(a1);
  }
  else {
    if (VERBOSE >1)
      fprintf(stderr," => CROSS %p/%p\n", a1, a2);
    accessBlock(a1);
    accessBlock(a2);
  }

  // consistency check
  if (DEBUG) {
    checkCount++;
    if (checkCount > 100000) {
      checkStack();
      checkCount = 0;
    }
  }
}

VOID memRead(ADDRINT addr, UINT32 size)
{
  if (VERBOSE >1)
    fprintf(stderr,"R %p/%d", (void*)addr, size);
  memAccess(addr, size);
}

VOID memWrite(ADDRINT addr, UINT32 size)
{
  if (VERBOSE >1)
    fprintf(stderr,"W %p/%d", (void*)addr, size);
  memAccess(addr, size);
}

VOID stackAccess()
{
  stackAccesses++;
}

/* ===================================================================== */
/* Instrumentation                                                       */
/* ===================================================================== */

VOID Instruction(INS ins, VOID* v)
{
  if (IGNORE_STACK && (INS_IsStackRead(ins) || INS_IsStackWrite(ins))) {
    INS_InsertPredicatedCall( ins, IPOINT_BEFORE, (AFUNPTR)stackAccess,
			      IARG_END);
    return;
  }

  UINT32 memOperands = INS_MemoryOperandCount(ins);
  for (UINT32 memOp = 0; memOp < memOperands; memOp++) {
    if (INS_MemoryOperandIsRead(ins, memOp))
      INS_InsertPredicatedCall( ins, IPOINT_BEFORE, (AFUNPTR)memRead,
				IARG_MEMORYOP_EA, memOp,
				IARG_MEMORYREAD_SIZE,
				IARG_END);

    if (INS_MemoryOperandIsWritten(ins, memOp))
      INS_InsertPredicatedCall( ins, IPOINT_BEFORE, (AFUNPTR)memWrite,
				IARG_MEMORYOP_EA, memOp,
				IARG_MEMORYREAD_SIZE,
				IARG_END);
  }
}


/* ===================================================================== */
/* Output results at exit                                                */
/* ===================================================================== */

VOID Exit(INT32 code, VOID *v)
{
  if (DEBUG)
    checkStack();

  unsigned long aCount = 0;
  for(const Bucket& b: buckets) {
    fprintf(stderr, "B [%d ...] %lu\n",
	    b.minDist(), b.aCount);
    aCount += b.aCount;
  }

  fprintf(stderr,
	  "Statistics:\n"
	  "  memory blocks accessed: %lu\n"
	  "  number of accesses:     %lu\n"
	  "  ignored stack accesses: %lu\n",
	  addrMap.size(), aCount, stackAccesses);
}

/* ===================================================================== */
/* Usage/Main Function of the Pin Tool                                   */
/* ===================================================================== */

INT32 Usage()
{
  PIN_ERROR( "PinDist: Get the Stack Reuse Distance Histogram\n" 
	     + KNOB_BASE::StringKnobSummary() + "\n");
  return -1;
}


int main (int argc, char *argv[])
{
  if (PIN_Init(argc, argv)) return Usage();

  // add buckets [0-1023], [1K - 2K-1], ... [1G - ]
  buckets.push_back( Bucket(0) );
  for(int i=1024; i< 128*1024*1024; i*=2)
    buckets.push_back( Bucket(i) );
  buckets.push_back( Bucket(0) ); // for "infinite" distance
  nextBucket = 1;
  stackAccesses = 0;

  PIN_InitSymbols();
  INS_AddInstrumentFunction(Instruction, 0);
  PIN_AddFiniFunction(Exit, 0);
  
  PIN_StartProgram();
  return 0;	
}
