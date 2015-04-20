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

#define DEBUG_STATS 0

typedef UINT32 u32;
typedef unsigned long long u64;
typedef unsigned long Addr;

#if DEBUG_STATS
#define MAX_PREDCOUNTERS 50
u64 predCounters[50];
#endif

struct AddrPredict {
  Addr currentAddr;
  Addr predictedAddr;
  INT32 stride;

  Addr iaddr;
  Addr lastAddr;
  INT32 lastStride;
  UINT32 missCount;
  UINT32 count;
};

#define MAXPREDS 1000000
struct AddrPredict apreds[MAXPREDS];
int nextAPred = 0;

#define MAX_EXITS 100000
int exit_from[MAX_EXITS];
int exit_to[MAX_EXITS];
u64 exit_counter[MAX_EXITS];
int nextExit = 0;

#define MAXDIFF 100000

#define likely(x) __builtin_expect((x),1)

void storeCurrentAddr(struct AddrPredict* p, Addr a)
{
  p->currentAddr = a;
}

__attribute__((always_inline)) inline
void check(struct AddrPredict* p, Addr a)
{
    if (likely(a == p->predictedAddr)) {
	p->predictedAddr += p->stride;

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
#if DEBUG_STATS
  predCounters[next-off]++;
#endif

  struct AddrPredict* p = &apreds[off];
  while(1) {
    check(p, p->currentAddr);
    off++;
    if (off == next) return;
    p++;
  }
}

template<int N>
void afterExitN(UINT32 off)
{
  struct AddrPredict* p;
  p = &apreds[off];
  check(p, p->currentAddr);
  for(int n=1; n<N; n++) {
    p++;
    check(p, p->currentAddr);
  }
}

void incExitCounter(int off)
{
  exit_counter[off]++;
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
    p->stride = 0;
    p->lastAddr = 0;
    p->lastStride = 0;
    p->iaddr = iaddr;
    p->missCount = 0;
    p->count = 0;

    nextAPred++;
    return nextAPred-1;
}

int getExit(int from, int to)
{
    if (nextExit >= MAX_EXITS) {
	fprintf(stderr, "ERROR: Too many trace exit counters needed\n");
	exit(1);
    }

    int c = nextExit;
    exit_from[c] = from;
    exit_to[c] = to;
    exit_counter[c] = 0;

    nextExit++;
    return c;
}

VOID Instruction(INS ins, VOID*)
{
    for (UINT32 i = 0; i < INS_MemoryOperandCount(ins); i++) {
	if (INS_MemoryOperandIsRead(ins,i) ||
	    INS_MemoryOperandIsWritten(ins,i)) {
	    int off = getAPred(INS_Address(ins));

	    INS_InsertCall( ins, IPOINT_BEFORE,
			    (AFUNPTR) storeCurrentAddr,
			    IARG_PTR, &apreds[off],
			    IARG_MEMORYOP_EA, i,
			    IARG_END);
	}
    }
}

AFUNPTR afterExitPtr(int n)
{
  if (n==0) return 0;
  //return (AFUNPTR) afterExit;

#if DEBUG_STATS
  return (AFUNPTR) afterExit;
#endif

  switch(n) {
  case 1: return (AFUNPTR) afterExitN<1>;
  case 2: return (AFUNPTR) afterExitN<2>;
  case 3: return (AFUNPTR) afterExitN<3>;
  case 4: return (AFUNPTR) afterExitN<4>;
  case 5: return (AFUNPTR) afterExitN<5>;
  case 6: return (AFUNPTR) afterExitN<6>;
  case 7: return (AFUNPTR) afterExitN<7>;
  case 8: return (AFUNPTR) afterExitN<8>;
  case 9: return (AFUNPTR) afterExitN<9>;
  }
  return (AFUNPTR) afterExit;
}

VOID Trace(TRACE trace, VOID *v)
{
    int first = nextAPred;

    BBL bbl;
    for (bbl = TRACE_BblHead(trace); BBL_Valid(bbl); bbl = BBL_Next(bbl)) {

        INS ins;
        for(ins= BBL_InsHead(bbl); INS_Valid(ins); ins = INS_Next(ins)) {
            Instruction(ins, 0);
        }

        if (INS_IsBranchOrCall(BBL_InsTail(bbl))) {
          AFUNPTR p = afterExitPtr(nextAPred-first);
	  if (p == (AFUNPTR) afterExit) {
	    BBL_InsertCall( bbl, IPOINT_TAKEN_BRANCH, p,
                            IARG_UINT32, first,
                            IARG_UINT32, nextAPred,
                            IARG_END);
	  }
	  else if (p != 0) {
	    BBL_InsertCall( bbl, IPOINT_TAKEN_BRANCH, p,
	                    IARG_UINT32, first,
                            IARG_END);
          }
	  BBL_InsertCall( bbl, IPOINT_TAKEN_BRANCH,
                          (AFUNPTR) incExitCounter,
                          IARG_UINT32, getExit(first, nextAPred),
                          IARG_END);
	}
    }
    if (TRACE_HasFallThrough(trace) && (nextAPred > first)) {
      AFUNPTR p = afterExitPtr(nextAPred-first);
      if (p == (AFUNPTR) afterExit) {
        TRACE_InsertCall( trace, IPOINT_AFTER, p,
			  IARG_UINT32, first,
			  IARG_UINT32, nextAPred,
                          IARG_END);
      }
      else if (p != 0) {
        TRACE_InsertCall( trace, IPOINT_AFTER, p,
			  IARG_UINT32, first,
                          IARG_END);
      }
      TRACE_InsertCall( trace, IPOINT_AFTER, (AFUNPTR) incExitCounter,
                        IARG_UINT32, getExit(first, nextAPred),
                        IARG_END);
    }
}

VOID Fini(int code, VOID * v)
{
    u64 totalMissCount = 0;
    u64 totalHitCount = 0;
    u64 totalCount = 0;
    u64 exitCounter = 0;

    for(int i=0; i<nextExit; i++) {
      for(int j=exit_from[i]; j<exit_to[i]; j++) {
	apreds[j].count += exit_counter[i];
	exitCounter += exit_counter[i];
      }
    }

    for(int i=0; i<nextAPred; i++) {
	struct AddrPredict* p = &apreds[i];

//	if (p1->hitCount + p2->missCount > 1000)
//	    fprintf(stderr, "%5i: I %lx, %d/%d\n",
//		    i, p2->iaddr, p1->hitCount, p1->hitCount + p2->missCount);

	totalMissCount += p->missCount;
	totalCount += p->count;
    }
    totalHitCount = totalCount - totalMissCount;

#if DEBUG_STATS
    for(int i=0; i<MAX_PREDCOUNTERS; i++) {
      if (predCounters[i] == 0) continue;
      fprintf(stderr, " PredCounters %2d : %llu\n", i, predCounters[i]);
    }
#endif

    fprintf(stderr, "Predicted: %llu / %llu = %5.2f %%\n",
            totalHitCount, totalCount,
	    100.0 * (double)totalHitCount / (double)totalCount);
    fprintf(stderr, " predictors: %d, trace exits: %d, exit counter: %llu\n",
            nextAPred, nextExit, exitCounter);
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

#if DEBUG_STATS
    for(int i=0; i<MAX_PREDCOUNTERS; i++)
      predCounters[i] = 0;
#endif

    // Never returns
    PIN_StartProgram();
    
    return 0;
}
