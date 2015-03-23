/* Pin Tool for writing memory access stream to file,
 * for later calculation of the Stack Reuse Distance Histogram
 *
 * (C) 2015, Josef Weidendorfer / LRR-TUM
 * GPLv2+ (see COPYING)
 */


#include "pin.H"
#include <stdio.h>
#include <cassert>

typedef void* Addr;

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

unsigned long stackAccesses, bCount, aCount;
FILE * trace;

/* ===================================================================== */
/* Handle Memory block access (aligned at multiple of MEMBLOCKLEN)       */
/* ===================================================================== */

void writeChunk(Addr* buf, int len)
{
#if MERGE_CHUNK
  sort(buf, buf+len);
#endif
  
  fwrite(buf, sizeof(Addr), len, trace);
  bCount += len;
}

void accessBlock(Addr a)
{
  static Addr mergeBuffer[CHUNKSIZE];
  static int ptr = 0;

  if (a == 0) {
    writeChunk(mergeBuffer, ptr);
    ptr = 0;
    return;
  }

  if (ptr < CHUNKSIZE) {
    mergeBuffer[ptr++] = a;
    return;
  }

  writeChunk(mergeBuffer, CHUNKSIZE);
  ptr = 0;
}


/* ===================================================================== */
/* Direct Callbacks                                                      */
/* ===================================================================== */

void memAccess(ADDRINT addr, UINT32 size)
{
  aCount++;

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
  fclose(trace);
  fprintf(stderr,
	  "Statistics:\n"
	  "  number of accesses:       %lu\n"
	  "  number of block accesses: %lu\n"
	  "  ignored stack accesses:   %lu\n",
	  aCount, bCount, stackAccesses);
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

  trace = fopen("trace.out", "w");
  stackAccesses = 0;
  aCount = 0;
  bCount = 0;

  PIN_InitSymbols();
  INS_AddInstrumentFunction(Instruction, 0);
  PIN_AddFiniFunction(Exit, 0);
  
  PIN_StartProgram();
  return 0;	
}
