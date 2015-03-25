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
#include <cmath>
#include <unistd.h>

#include "dist.cpp"

// Consistency checks?
#define DEBUG 0

// 2: Huge amount of debug output, 1: checks, 0: silent
#define VERBOSE 0

// uses INS_IsStackRead/Write: misleading with -fomit-frame-pointer
#define IGNORE_STACK 1

// collect addresses in chunk buffer before? (always worse)
#define MERGE_CHUNK 0
#define CHUNKSIZE 4096

// must be a power-of-two
#define MEMBLOCKLEN 64

unsigned long stackAccesses;

/* ===================================================================== */
/* Command line options                                                  */
/* ===================================================================== */

KNOB<int> KnobMinDist(KNOB_MODE_WRITEONCE, "pintool",
    "m", "4096", "minimum bucket distance");
KNOB<int> KnobDoubleSteps(KNOB_MODE_WRITEONCE, "pintool",
    "s", "1", "number of buckets for doubling distance");
KNOB<bool> KnobPIDPrefix(KNOB_MODE_WRITEONCE, "pintool",
    "p", "0", "prepend output by --PID--");

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

char* format(unsigned long aCount, unsigned long maxCount)
{
  static char out[20];

  if (maxCount > 999999999)
    sprintf(out, "%7.3f G", (double) aCount / 1000000000.0);
  else if (maxCount > 999999)
    sprintf(out, "%7.3f M", (double) aCount / 1000000.0);
  else if (aCount > 999)
    sprintf(out, "%3d %3d  ",
	    (int) aCount / 1000, (int)(aCount % 1000));
  else
    sprintf(out, "    %3d  ", (int) aCount);

  return out;
}

VOID Exit(INT32 code, VOID *v)
{
  int b, bNext;
  unsigned int min;
  unsigned long aCount, maxCount;
  unsigned long stack_size;
  char bar[] = "##############################################################";
  int barCount = strlen(bar), barSize;
  char pStr[20], bStr[20];

  if (KnobPIDPrefix.Value())
    sprintf(pStr, "--%5d-- ", getpid());
  else
    pStr[0] = 0;

  fprintf(stderr, "%sHistogram:\n", pStr);

  maxCount = 0;
  b = 1;
  do {
    bNext = RD_get_hist(b, min, aCount);
    if (aCount > maxCount) maxCount = aCount;
    b = bNext;
  } while(b!=0);

  bNext = RD_get_hist(0, min, aCount);
  fprintf(stderr, "%s[%8.3f MB ..] %s ==>\n",
	  pStr, (double)(min * MEMBLOCKLEN)/1000000.0,
	  format(aCount,aCount));
  b = bNext;
  do {
    bNext = RD_get_hist(b, min, aCount);
    barSize = (int)((double)aCount/maxCount*barCount);
    if (barSize > barCount) barSize = barCount;
    if (min>0)
      sprintf(bStr, "%8.3f MB ..",
	      (double)(min * MEMBLOCKLEN) / 1024.0 / 1024.0);
    else
      sprintf(bStr, "   inf/cold   ");
    fprintf(stderr, "%s[%s] %s %s\n",
	    pStr, bStr,
	    format(aCount, maxCount),
	    bar + barCount - barSize);
    b = bNext;
  } while(b!=0);


  RD_stat(stack_size, aCount);

  fprintf(stderr,
	  "%sStatistics:\n"
	  "%s  memory blocks accessed: %lu (%3.3f MB)\n"
	  "%s  number of accesses:     %lu\n"
	  "%s  ignored stack accesses: %lu\n",
	  pStr, pStr, stack_size, ((double)stack_size * MEMBLOCKLEN)/1000000.0,
	  pStr, aCount,
	  pStr, stackAccesses);
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
