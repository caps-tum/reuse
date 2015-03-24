/* Pin Tool for
 * calculation of the Stack Reuse Distance Histogram
 *
 * (C) 2015, Josef Weidendorfer / LRR-TUM
 * GPLv2+ (see COPYING)
 */


#include "pin.H"

#include <stdio.h>
#include <cassert>
#include <cstring>
#include <unistd.h>

#include "dist.cpp"

// Consistency checks?
#define DEBUG 0

// 2: Huge amount of debug output, 1: checks, 0: silent
#define VERBOSE 0

// uses INS_IsStackRead/Write: misleading with -fomit-frame-pointer
#define IGNORE_STACK 1

// collect addresses in chunk buffer before?
#define MERGE_CHUNK 0
#define CHUNKSIZE 4096

// must be a power-of-two
#define MEMBLOCKLEN 64

unsigned long stackAccesses;

/* ===================================================================== */
/* Handle Memory block access (aligned at multiple of MEMBLOCKLEN)       */
/* ===================================================================== */

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
    RD_accessBlock(mergeBuffer[ptr]);
  }
  ptr = 0;
}
#define RD_accessBlock accessMerging
#endif


/* ===================================================================== */
/* Direct Callbacks                                                      */
/* ===================================================================== */

void memAccess(ADDRINT addr, UINT32 size)
{
  Addr a1 = (void*) (addr & ~(MEMBLOCKLEN-1));
  Addr a2 = (void*) ((addr+size-1) & ~(MEMBLOCKLEN-1));
  if (a1 == a2) {
    if (VERBOSE >1)
      fprintf(stderr," => %p\n", a1);
    RD_accessBlock(a1);
  }
  else {
    if (VERBOSE >1)
      fprintf(stderr," => CROSS %p/%p\n", a1, a2);
    RD_accessBlock(a1);
    RD_accessBlock(a2);
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
        IARG_UINT32, INS_MemoryOperandSize(ins, memOp),
				IARG_END);

    if (INS_MemoryOperandIsWritten(ins, memOp))
      INS_InsertPredicatedCall( ins, IPOINT_BEFORE, (AFUNPTR)memWrite,
				IARG_MEMORYOP_EA, memOp,
        IARG_UINT32, INS_MemoryOperandSize(ins, memOp),
				IARG_END);
  }
}


/* ===================================================================== */
/* Output results at exit                                                */
/* ===================================================================== */

VOID Exit(INT32 code, VOID *v)
{
  int b, bNext;
  unsigned int min;
  unsigned long aCount, maxCount;
  unsigned long stack_size;
  char bar[] = "##############################################################";
  int barCount = strlen(bar), barSize;

  fprintf(stderr, "Histogram (PID %d):\n", getpid());

  maxCount = 0;
  b = 1;
  do {
    bNext = RD_get_hist(b, min, aCount);
    if (aCount > maxCount) maxCount = aCount;
    b = bNext;
  } while(b!=0);

  bNext = RD_get_hist(0, min, aCount);
  fprintf(stderr, "[%8.3f MB .. ] %8lu ==>\n",
      (double)(min * MEMBLOCKLEN)/1000000.0, aCount);
  b = bNext;
  do {
    bNext = RD_get_hist(b, min, aCount);
    barSize = (int)((double)aCount/maxCount*barCount);
    if (barSize > barCount) barSize = barCount;
    fprintf(stderr, "[%8.3f MB .. ] %8lu %s\n",
      (double)(min * MEMBLOCKLEN)/1000000.0,
      aCount,
      bar + barCount - barSize);
    b = bNext;
  } while(b!=0);


  RD_stat(stack_size, aCount);

  fprintf(stderr,
	  "Statistics:\n"
	  "  memory blocks accessed: %lu (%3.3f MB)\n"
	  "  number of accesses:     %lu\n"
	  "  ignored stack accesses: %lu\n",
	  stack_size, ((double)stack_size * MEMBLOCKLEN)/1000000.0,
	  aCount, stackAccesses);
}

/* ===================================================================== */
/* Usage/Main Function of the Pin Tool                                   */
/* ===================================================================== */

KNOB<int> KnobMinDist(KNOB_MODE_WRITEONCE, "pintool",
    "m", "4096", "minimum bucket distance");
KNOB<int> KnobDoubleSteps(KNOB_MODE_WRITEONCE, "pintool",
    "s", "1", "number of buckets for doubling distance");


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
  double d = KnobMinDist.Value();
  int s = KnobDoubleSteps.Value();
  double f = pow(2, 1.0/s);
  RD_init((int)(d / MEMBLOCKLEN));
  for(d*=f; d< 1024*1024*1024; d*=f)
    RD_addBucket((int)(d / MEMBLOCKLEN));

  stackAccesses = 0;

  PIN_InitSymbols();
  INS_AddInstrumentFunction(Instruction, 0);
  PIN_AddFiniFunction(Exit, 0);
  
  PIN_StartProgram();
  return 0;	
}
