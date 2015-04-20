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
#include <vector>

typedef UINT32 u32;
typedef unsigned long long u64;
typedef unsigned long Addr;

struct AddrPredict {
  Addr currentAddr;
  Addr predictedAddr;
  UINT32 hitCount;
  INT32 stride;

  Addr iaddr;
  Addr lastAddr;
  INT32 lastStride;
  UINT32 missCount;
};

#define MAXPREDS 1000000
struct AddrPredict apreds[MAXPREDS];
int nextAPred = 0;

#define MAXDIFF 100000

#define likely(x) __builtin_expect((x),1)

void storeCurrentAddr(struct AddrPredict* p, Addr a)
{
  p->currentAddr = a;
}

void check(struct AddrPredict* p, Addr a)
{
    if (likely(a == p->predictedAddr)) {
	p->predictedAddr += p->stride;
	p->hitCount++;
//	fprintf(stderr, "H I %lx %4d/%-4d: %lx      => %lx/%d\n",
//		p->iaddr, p->hitCount, p->missCount + p->hitCount,
//		a, p->addr, p->stride);
	return;
    }

    p->missCount++;

    long diff = a - p->lastAddr;
    p->lastAddr = a;

    if (p->lastStride == diff) {
	p->stride = diff;
	p->predictedAddr = a + diff;
    }
    else {
	p->lastStride = (diff > -MAXDIFF && diff < MAXDIFF) ? diff : 0;
	p->predictedAddr = a + p->stride;
    }
    //	if (p->missCount > 1000)
    //	    fprintf(stderr, "M I %lx %4d/%-4d: %lx %4d => %lx/%d\n",
    //		    p->iaddr, p->hitCount, p->missCount + p->hitCount,
    //		    p->lastAddr, p->lastStride,
    //		    p->addr, p->stride);
}

void afterExit(UINT32 off, UINT32 next)
{
  while(off < next) {
    struct AddrPredict* p = &apreds[off];
    check(p, p->currentAddr);
    off++;
  }
}

/* ===================================================================== */
// Callbacks

int getAPred(Addr iaddr)
{
    if (nextAPred >= MAXPREDS) {
	fprintf(stderr, "ERROR: Too many predictors needed\n");
	exit(1);
    }

    struct AddrPredict* p = &apreds[nextAPred];

    p->predictedAddr = 0;
    p->hitCount = 0;
    p->stride = 0;
    p->lastAddr = 0;
    p->lastStride = 0;
    p->iaddr = iaddr;
    p->missCount = 0;

    nextAPred++;
    return nextAPred-1;
}


VOID Instruction(INS ins, VOID*)
{
    for (UINT32 i = 0; i < INS_MemoryOperandCount(ins); i++) {
	if (INS_MemoryOperandIsRead(ins,i) || INS_MemoryOperandIsWritten(ins,i)) {
	    int off = getAPred(INS_Address(ins));

	    INS_InsertCall( ins, IPOINT_BEFORE,
			    (AFUNPTR) storeCurrentAddr,
			    IARG_PTR, &apreds[off],
			    IARG_MEMORYOP_EA, i,
			    IARG_END);
	}
    }
}


VOID Trace(TRACE trace, VOID *v)
{
    int firstAPred = nextAPred;

    BBL bbl;
    for (bbl = TRACE_BblHead(trace); BBL_Valid(bbl); bbl = BBL_Next(bbl)) {

        INS ins;
        for(ins= BBL_InsHead(bbl); INS_Valid(ins); ins = INS_Next(ins)) {
            Instruction(ins, 0);
        }

        if (INS_IsBranchOrCall(BBL_InsTail(bbl)) && (nextAPred > firstAPred))
	    BBL_InsertCall( bbl, IPOINT_TAKEN_BRANCH, (AFUNPTR) afterExit,
                            IARG_UINT32, firstAPred,
                            IARG_UINT32, nextAPred,
                            IARG_END);

    }
    if (TRACE_HasFallThrough(trace) && (nextAPred > firstAPred))
        TRACE_InsertCall( trace, IPOINT_AFTER, (AFUNPTR) afterExit,
			  IARG_UINT32, firstAPred,
			  IARG_UINT32, nextAPred,
                          IARG_END);
}

VOID Fini(int code, VOID * v)
{
    u64 totalHitCount = 0;
    u64 totalCount = 0;

    for(int i=0; i<nextAPred; i++) {
	struct AddrPredict* p = &apreds[i];

//	if (p1->hitCount + p2->missCount > 1000)
//	    fprintf(stderr, "%5i: I %lx, %d/%d\n",
//		    i, p2->iaddr, p1->hitCount, p1->hitCount + p2->missCount);

	totalHitCount += p->hitCount;
	totalCount += p->missCount + p->hitCount;
    }

    fprintf(stderr, "Predicted: %llu / %llu (preds: %d) = %5.2f %%\n",
	    totalHitCount, totalCount, nextAPred,
	    100.0 * (double)totalHitCount / (double)totalCount);
}

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

    TRACE_AddInstrumentFunction(Trace, 0);
    //INS_AddInstrumentFunction(Instruction, 0);
    PIN_AddFiniFunction(Fini, 0);

    //predReg = PIN_ClaimToolRegister();

    // Never returns
    PIN_StartProgram();
    
    return 0;
}
