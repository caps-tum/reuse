/*
 * Dynamic Prediction for Memory Accesses
 */

#include "pin.H"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <map>


typedef UINT32 u32;
typedef unsigned long long u64;
typedef unsigned long Addr;

void doRead(THREADID t, Addr addr, UINT32 size)
{
  fprintf(stderr, "R T%2d A%lx S%d\n", t, addr, size);
}

void doWrite(THREADID t, Addr addr, UINT32 size)
{
  fprintf(stderr, "W T%2d A%lx S%d\n", t, addr, size);
}

/* ===================================================================== */
// Callbacks

VOID Instruction(INS ins, VOID*)
{
  for (UINT32 i = 0; i < INS_MemoryOperandCount(ins); i++) {
    int dSize = INS_MemoryOperandSize(ins, i);
    
    // TODO: Predicated moves
    if (INS_MemoryOperandIsRead(ins,i)) {
      INS_InsertCall( ins, IPOINT_BEFORE,
              (AFUNPTR) doRead,
		      IARG_THREAD_ID,
		      IARG_MEMORYOP_EA, i,
		      IARG_UINT32, dSize,
		      IARG_END);
    }
    else if (INS_MemoryOperandIsWritten(ins,i)) {
      INS_InsertCall( ins, IPOINT_BEFORE,
              (AFUNPTR) doWrite,
		      IARG_THREAD_ID,
		      IARG_MEMORYOP_EA, i,
		      IARG_UINT32, dSize,
		      IARG_END);
    }
  }
}

VOID Fini(int code, VOID * v)
{}

/* ===================================================================== */
/* Main                                                                  */
/* ===================================================================== */

int main(int argc, char *argv[])
{
    PIN_InitSymbols();
    if( PIN_Init(argc,argv) ) {
      fprintf(stderr, "PinPredict V0.1\n");
      return 0;
    }

    INS_AddInstrumentFunction(Instruction, 0);
    PIN_AddFiniFunction(Fini, 0);

    // Never returns
    PIN_StartProgram();
    
    return 0;
}
