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

#define ONE_PRED 1

#if ONE_PRED
struct AddrPredict1 {
  Addr addr;
  UINT32 hitCount;
  INT32 stride;

  Addr iaddr;
  Addr lastAddr;
  INT32 lastStride;
  UINT32 missCount;
};
#else
struct AddrPredict1 {
  Addr addr;
  INT32 stride;
  UINT32 hitCount;
};

struct AddrPredict2 {
  Addr iaddr;
  Addr lastAddr;
  INT32 lastStride;
  UINT32 missCount;
};
#endif

#if 0
vector<struct AddrPredict> predictVector;
u64 p1Count, p2Count;

REG predReg;

void doRead(THREADID t, Addr addr, UINT32 size)
{
  fprintf(stderr, "R T%2d A%lx S%d\n", t, addr, size);
}

void doWrite(THREADID t, Addr addr, UINT32 size)
{
  fprintf(stderr, "W T%2d A%lx S%d\n", t, addr, size);
}

int checkPrediction(struct AddrPredict* p, Addr a)
{
  return (a == p->addr);
}

int isTrue(ADDRINT v)
{
  return v;
}

int isFalse(ADDRINT v)
{
  return !v;
}

int checkNotPrediction(struct AddrPredict* p, Addr a)
{
  return (a != p->addr);
}

void predictionOK(struct AddrPredict* p, Addr a)
{
  p1Count++;
}

void predictionNotOK(struct AddrPredict* p, Addr a)
{
  p2Count++;
}
#endif

#define MAXPREDS 1000000
#if ONE_PRED
struct AddrPredict1 apreds1[MAXPREDS];
#else
struct AddrPredict1 apreds1[MAXPREDS];
struct AddrPredict2 apreds2[MAXPREDS];
#endif
int nextAPred = 0;

#define MAXDIFF 100000

#define likely(x) __builtin_expect((x),1)

void check(struct AddrPredict1* p1, Addr a)
{
    if (likely(a == p1->addr)) {
	//if (p1->stride)
	p1->addr += p1->stride;
	p1->hitCount++;
//	fprintf(stderr, "H I %lx %4d/%-4d: %lx      => %lx/%d\n",
//		p->iaddr, p->hitCount, p->missCount + p->hitCount,
//		a, p->addr, p->stride);
	return;
    }

#if ONE_PRED
    struct AddrPredict1* p2 = p1;
#else
    int off = p1 - apreds1;
    struct AddrPredict2* p2 = apreds2 + off;
//    struct AddrPredict2* p2 = (struct AddrPredict2*) (((char*)p1 -
//						       (char*)apreds1) +
//						      (char*)apreds2);
#endif
    p2->missCount++;

    long diff = a - p2->lastAddr;
    p2->lastAddr = a;

    if (p2->lastStride == diff) {
	p1->stride = diff;
	p1->addr = a + diff;
    }
    else {
	p2->lastStride = (diff > -MAXDIFF && diff < MAXDIFF) ? diff : 0;
	p1->addr = a + p1->stride;
    }
    //	if (p->missCount > 1000)
    //	    fprintf(stderr, "M I %lx %4d/%-4d: %lx %4d => %lx/%d\n",
    //		    p->iaddr, p->hitCount, p->missCount + p->hitCount,
    //		    p->lastAddr, p->lastStride,
    //		    p->addr, p->stride);
}


/* ===================================================================== */
// Callbacks

int getAPred(Addr iaddr)
{
    if (nextAPred >= MAXPREDS) {
	fprintf(stderr, "ERROR: Too many predictors needed\n");
	exit(1);
    }

    struct AddrPredict1* p1 = &apreds1[nextAPred];
#if ONE_PRED
    struct AddrPredict1* p2 = p1;
#else
    struct AddrPredict2* p2 = &apreds2[nextAPred];
#endif


    p1->addr = 0;
    p1->hitCount = 0;
    p1->stride = 0;

    p2->lastAddr = 0;
    p2->lastStride = 0;
    p2->iaddr = iaddr;
    p2->missCount = 0;

//    fprintf(stderr, "  I %lx\n", p->iaddr);

    nextAPred++;
    return nextAPred-1;
}


VOID Instruction(INS ins, VOID*)
{
  for (UINT32 i = 0; i < INS_MemoryOperandCount(ins); i++) {
    //int dSize = INS_MemoryOperandSize(ins, i);
    
    if (INS_MemoryOperandIsRead(ins,i) || INS_MemoryOperandIsWritten(ins,i)) {
	int off = getAPred(INS_Address(ins));
#if 1
      INS_InsertCall( ins, IPOINT_BEFORE,
	(AFUNPTR) check,
	IARG_PTR, &apreds1[off],
	IARG_MEMORYOP_EA, i,
	IARG_END);
#endif
#if 0
      INS_InsertCall( ins, IPOINT_BEFORE,
		      (AFUNPTR) checkPrediction,
		      IARG_PTR, &p,
		      IARG_MEMORYOP_EA, i,
		      IARG_RETURN_REGS, predReg,
		      IARG_END);
      INS_InsertIfCall( ins, IPOINT_BEFORE,
			(AFUNPTR) isTrue,
			IARG_REG_VALUE, predReg,
			IARG_END);
      INS_InsertThenCall( ins, IPOINT_BEFORE,
			  (AFUNPTR) predictionOK,
			  IARG_PTR, &p,
			  IARG_MEMORYOP_EA, i,
			  IARG_END);
      INS_InsertIfCall( ins, IPOINT_BEFORE,
			(AFUNPTR) isFalse,
			IARG_REG_VALUE, predReg,
			IARG_END);
      INS_InsertThenCall( ins, IPOINT_BEFORE,
			  (AFUNPTR) predictionNotOK,
			  IARG_PTR, &p,
			  IARG_MEMORYOP_EA, i,
			  IARG_END);
#endif
#if 0
      INS_InsertIfCall( ins, IPOINT_BEFORE,
			(AFUNPTR) checkPrediction,
			IARG_PTR, &p,
			IARG_MEMORYOP_EA, i,
			IARG_END);
      INS_InsertThenCall( ins, IPOINT_BEFORE,
			  (AFUNPTR) predictionOK,
			  IARG_PTR, &p,
			  IARG_MEMORYOP_EA, i,
			  IARG_END);
      INS_InsertIfCall( ins, IPOINT_BEFORE,
			(AFUNPTR) checkNotPrediction,
			IARG_PTR, &p,
			IARG_MEMORYOP_EA, i,
			IARG_END);
      INS_InsertThenCall( ins, IPOINT_BEFORE,
			  (AFUNPTR) predictionNotOK,
			  IARG_PTR, &p,
			  IARG_MEMORYOP_EA, i,
			  IARG_END);
#endif
    }
  }
}

VOID Fini(int code, VOID * v)
{
    u64 totalHitCount = 0;
    u64 totalCount = 0;

    for(int i=0; i<nextAPred; i++) {
	struct AddrPredict1* p1 = &apreds1[i];
#if ONE_PRED
	struct AddrPredict1* p2 = p1;
#else
	struct AddrPredict2* p2 = &apreds2[i];
#endif
//	if (p1->hitCount + p2->missCount > 1000)
//	    fprintf(stderr, "%5i: I %lx, %d/%d\n",
//		    i, p2->iaddr, p1->hitCount, p1->hitCount + p2->missCount);

	totalHitCount += p1->hitCount;
	totalCount += p2->missCount + p1->hitCount;
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

    INS_AddInstrumentFunction(Instruction, 0);
    PIN_AddFiniFunction(Fini, 0);

    //predReg = PIN_ClaimToolRegister();

    // Never returns
    PIN_StartProgram();
    
    return 0;
}
