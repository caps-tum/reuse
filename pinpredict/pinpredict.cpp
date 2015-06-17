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
#include <algorithm>

#define DEBUG_STATS 0

#define IGNORE_STACK 1

typedef UINT32 u32;
typedef unsigned long long u64;
typedef unsigned long Addr;

struct AddrPredict {
  Addr currentAddr;
  Addr predictedAddr;
  INT32 predictedStride;

  INT32 lastStride;
  Addr lastAddr;
  UINT32 missCount;
  UINT32 count;
};

struct PredInfo {
  Addr iaddr;
  string funcname;
  string filename;
  int line, column;
};

#define MAX_PREDS 100000
struct PredInfo predInfo[MAX_PREDS];
int nextAPred = 0;

#define MAX_EXITS 50000
int exit_from[MAX_EXITS];
int exit_to[MAX_EXITS];
int exit_stackAccs[MAX_EXITS];
int nextExit = 0;

// for DEBUG_STATS
#define MAX_PREDCOUNTERS 50

struct ThreadData {
    struct AddrPredict apreds[MAX_PREDS];
    u64 exit_counter[MAX_EXITS];
    u64 predCounters[MAX_PREDCOUNTERS];
};

#define MAX_THREADS 50
struct ThreadData* tdata[MAX_THREADS];

REG predReg;

#define MAXDIFF 100000

#define likely(x)   __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)

/* ===================================================================== */
// Callbacks

void* returnPredPtr(THREADID tid, UINT off)
{
    return &(tdata[tid]->apreds[off]);
}

template<int N>
void* returnPredPtrN(struct AddrPredict* p)
{
    return p+N;
}

template<int N>
void storeAddrN(struct AddrPredict* p, Addr a)
{
    (p+N)->currentAddr = a;
}


__attribute__((always_inline)) inline
void check(struct AddrPredict* p, Addr a)
{
    if (likely(a == p->predictedAddr)) {
#if 0
	// not resetting lastStride to predicted on successful prediction
	// can sometimes result in wrong prediction
	// however, it does not change results much and reduces performance
	if (unlikely(p->predictedStride != p->lastStride))
	    p->lastStride = p->predictedStride;
#endif
	p->predictedAddr += p->predictedStride;

//	fprintf(stderr, "H I %lx %4d/%-4d: %lx      => %lx/%d\n",
//		p->iaddr, p->hitCount, p->missCount + p->hitCount,
//		a, p->addr, p->stride);
	return;
    }

    p->missCount++;

    long stride = a - p->lastAddr;
    p->lastAddr = a;

    if (p->lastStride == stride) {
	p->predictedStride = stride;
	p->predictedAddr = a + stride;
    }
    else {
	p->lastStride = (stride > -MAXDIFF && stride < MAXDIFF) ? stride : 0;
	p->predictedAddr = a + p->predictedStride;
    }
    //	if (p->missCount > 1000)
    //	    fprintf(stderr, "M I %lx %4d/%-4d: %lx %4d => %lx/%d\n",
    //		    p->iaddr, p->hitCount, p->missCount + p->hitCount,
    //		    p->lastAddr, p->lastStride,
    //		    p->addr, p->stride);
}

void afterExit(THREADID tid, UINT32 off, UINT32 next)
{
    //fprintf(stderr, "afterExit tid %d, off %d, next %d\n", tid, off, next);

#if DEBUG_STATS
  tdata[tid]->predCounters[next-off]++;
#endif
  if (off == next) return;

  struct AddrPredict* p = &(tdata[tid]->apreds[off]);
  while(1) {
    check(p, p->currentAddr);
    off++;
    if (off == next) return;
    p++;
  }
}

template<int N>
void afterExitN(struct AddrPredict* p)
{
  check(p, p->currentAddr);
  for(int n=1; n<N; n++) {
    p++;
    check(p, p->currentAddr);
  }
}

void incExitCounter(THREADID tid, int off)
{
  tdata[tid]->exit_counter[off]++;
}


/* ===================================================================== */
// low frequency callbacks: thread creation etc.

VOID ThreadStart(THREADID tid, CONTEXT *ctxt, INT32 flags, VOID *v)
{
    if (tid >= MAX_THREADS) {
      fprintf(stderr, "ERROR: too many threads (tid %d)\n", tid);
      exit(1);
    }
    //fprintf(stderr, "== ThreadStart %d\n", tid);

    struct ThreadData* d;
    d = (struct ThreadData*) malloc(sizeof(struct ThreadData));

    for(int i=0; i<MAX_PREDS; i++) {
        struct AddrPredict* p = &(d->apreds[i]);

        p->currentAddr = 0;
        p->predictedAddr = 0;
	p->predictedStride = 0;
        p->lastAddr = 0;
        p->lastStride = 0;
        p->missCount = 0;
        p->count = 0;
    }

    for(int i=0; i<MAX_EXITS; i++)
        d->exit_counter[i] = 0;

    for(int i=0; i<MAX_PREDCOUNTERS; i++)
        d->predCounters[i] = 0;

    tdata[tid] = d;
}


/* ===================================================================== */
// Instrumentation

int getAPred(Addr iaddr)
{
    if (nextAPred >= MAX_PREDS) {
	fprintf(stderr, "ERROR: Too many predictors needed\n");
	exit(1);
    }
    int p = nextAPred;
    predInfo[p].iaddr = iaddr;
    predInfo[p].funcname = RTN_FindNameByAddress(iaddr);
    PIN_GetSourceLocation(iaddr,
			  &(predInfo[p].column),
			  &(predInfo[p].line),
			  &(predInfo[p].filename));

    nextAPred++;
    return p;
}

int getExit(int from, int to, int stackAccs)
{
    if (nextExit >= MAX_EXITS) {
	fprintf(stderr, "ERROR: Too many trace exit counters needed\n");
	exit(1);
    }

    int c = nextExit;
    exit_from[c] = from;
    exit_to[c] = to;
    exit_stackAccs[c] = stackAccs;

    nextExit++;
    return c;
}

#define MAX_STOREFUNC 30
AFUNPTR storeFuncs[MAX_STOREFUNC];

template<int N>
void storeFuncGen()
{
    storeFuncs[N] = (AFUNPTR) storeAddrN<N>;
    storeFuncGen<N-1>();
}

template<>
void storeFuncGen<0>() { storeFuncs[0] = (AFUNPTR) storeAddrN<0>; }


AFUNPTR storeAddrPtr(int n)
{
    if (n>=0 && n<MAX_STOREFUNC) return storeFuncs[n];

    fprintf(stderr, "ERROR: in address store instrumentation\n");
    exit(1);
}

// up to which offset store into buffer for current address
#define MAX_ADDROFFSET 30

// for instrumentation: number of accesses already observed in current trace
int accNumber = 0;
int accBase = 0;
int accStack = 0;

VOID Instruction(INS ins, VOID*)
{
#if IGNORE_STACK
    // predictability of stack accesses not interesting: always from L1
    if (INS_IsStackRead(ins) || INS_IsStackWrite(ins)) {
        accStack++;
        return;
    }
#endif

    for (UINT32 i = 0; i < INS_MemoryOperandCount(ins); i++) {
        if (INS_MemoryOperandIsRead(ins,i) ||
                INS_MemoryOperandIsWritten(ins,i)) {
	  unsigned long int a = INS_Address(ins);
	  int off = getAPred(a);

            if (accNumber == 0)
                INS_InsertCall( ins, IPOINT_BEFORE, (AFUNPTR) returnPredPtr,
                                IARG_THREAD_ID,
                                IARG_UINT32, off,
                                IARG_RETURN_REGS, predReg,
                                IARG_END);

            if (accNumber - accBase >= MAX_ADDROFFSET) {
                INS_InsertCall( ins, IPOINT_BEFORE,
                                (AFUNPTR) returnPredPtrN<MAX_ADDROFFSET>,
                                IARG_REG_VALUE, predReg,
                                IARG_RETURN_REGS, predReg,
                                IARG_END);
                accBase += MAX_ADDROFFSET;
            }

            AFUNPTR p = storeAddrPtr(accNumber - accBase);
            INS_InsertCall( ins, IPOINT_BEFORE, p,
                            IARG_REG_VALUE, predReg,
                            IARG_MEMORYOP_EA, i,
                            IARG_END);

            accNumber++;
        }
    }
}

AFUNPTR afterExitPtr(int n)
{
#if DEBUG_STATS
  return (AFUNPTR) afterExit;
#endif

  if (n==0) return 0;

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
    accNumber = 0;
    accBase = 0;
    accStack = 0;

    BBL bbl;
    for (bbl = TRACE_BblHead(trace); BBL_Valid(bbl); bbl = BBL_Next(bbl)) {

        INS ins;
        for(ins= BBL_InsHead(bbl); INS_Valid(ins); ins = INS_Next(ins)) {
            Instruction(ins, 0);
        }

        if (INS_IsBranchOrCall(BBL_InsTail(bbl))) {
	    AFUNPTR p = afterExitPtr(nextAPred-first);
	    //AFUNPTR p = (AFUNPTR)afterExit;
            if (p == (AFUNPTR) afterExit) {
                BBL_InsertCall( bbl, IPOINT_TAKEN_BRANCH, p,
                                IARG_THREAD_ID,
                                IARG_UINT32, first,
                                IARG_UINT32, nextAPred,
                                IARG_END);
            }
            else if (p != 0) {
                // call to afterExitN
                if (accBase > 0) {
                    // predReg has to be set again
                    BBL_InsertCall( bbl, IPOINT_TAKEN_BRANCH,
                                    (AFUNPTR) returnPredPtr,
                                    IARG_THREAD_ID,
                                    IARG_UINT32, first,
                                    IARG_RETURN_REGS, predReg,
                                    IARG_END);
                }
                BBL_InsertCall( bbl, IPOINT_TAKEN_BRANCH, p,
                                IARG_REG_VALUE, predReg,
                                IARG_END);
            }
            BBL_InsertCall( bbl, IPOINT_TAKEN_BRANCH,
                            (AFUNPTR) incExitCounter,
                            IARG_THREAD_ID,
                            IARG_UINT32, getExit(first, nextAPred, accStack),
                            IARG_END);
        }
    }
    if (TRACE_HasFallThrough(trace)) {
        AFUNPTR p = afterExitPtr(nextAPred-first);
        if (p == (AFUNPTR) afterExit) {
            TRACE_InsertCall( trace, IPOINT_AFTER, p,
                              IARG_THREAD_ID,
                              IARG_UINT32, first,
                              IARG_UINT32, nextAPred,
                              IARG_END);
        }
        else if (p != 0) {
            // call to afterExitN
            if (accBase > 0) {
                // predReg has to be set again
                TRACE_InsertCall( trace, IPOINT_AFTER,
                                  (AFUNPTR) returnPredPtr,
                                  IARG_THREAD_ID,
                                  IARG_UINT32, first,
                                  IARG_RETURN_REGS, predReg,
                                  IARG_END);
            }
            TRACE_InsertCall( trace, IPOINT_AFTER, p,
                              IARG_REG_VALUE, predReg,
                              IARG_END);
        }
        TRACE_InsertCall( trace, IPOINT_AFTER, (AFUNPTR) incExitCounter,
                          IARG_THREAD_ID,
                          IARG_UINT32, getExit(first, nextAPred, accStack),
                          IARG_END);
    }
    //fprintf(stderr, "Processed Trace %lx\n", TRACE_Address(trace));
}

struct namesorter {
    namesorter(const map<string,u64>& m):_m(m) {}
    bool operator()(string i, string j) { return _m.at(i)<_m.at(j);}
    const map<string,u64>& _m;
};

VOID Fini(int code, VOID * v)
{
    u64 totalMissCount = 0;
    u64 totalHitCount = 0;
    u64 totalCount = 0;
    u64 totalCount2 = 0;
    int accs;

    map<string,u64> routineTotals;
    map<string,u64> routineMisses;
    map<u64,string> sortedNames;

    fprintf(stderr,
            "==\n"
            "PinTool Results\n"
            "  predictors: %d, exits: %d\n",
            nextAPred, nextExit);

    for(int tid=0; tid<MAX_THREADS; tid++) {
        struct ThreadData* d = tdata[tid];
        if (d == 0) continue;

        u64 threadMissCount = 0;
        u64 threadHitCount = 0;
        u64 threadCount = 0;
        u64 threadCount2 = 0;
        int threadExitsUsed = 0;

        for(int i=0; i<nextExit; i++) {
            u64 counter = d->exit_counter[i];
            if (counter == 0) continue;
            accs = exit_to[i] - exit_from[i] + exit_stackAccs[i];
            threadCount2 += accs * counter;
            threadExitsUsed++;

            for(int j=exit_from[i]; j<exit_to[i]; j++) {
                d->apreds[j].count += counter;
                threadCount2 += counter;
            }
        }

	routineMisses.clear();
	routineTotals.clear();
	sortedNames.clear();
        for(int i=0; i<nextAPred; i++) {
            struct AddrPredict* p = &(d->apreds[i]);

            threadMissCount += p->missCount;
            threadCount += p->count;

	    routineTotals[predInfo[i].funcname] += p->count;
	    routineMisses[predInfo[i].funcname] += p->missCount;
        }
        threadHitCount = threadCount - threadMissCount;

        totalMissCount += threadMissCount;
        totalHitCount += threadHitCount;
        totalCount += threadCount;
        totalCount2 += threadCount2;


        fprintf(stderr, "  Thread %d\n", tid);
        fprintf(stderr, "    predicted%s: %llu / %llu = %5.2f %%\n",
                (threadCount == threadCount2) ? "":" (no stack accesses)",
                threadHitCount, threadCount,
                100.0 * (double)threadHitCount / (double)threadCount);
        fprintf(stderr, "    (exits used: %d, all accesses: %llu)\n",
                threadExitsUsed, threadCount2);

	vector<string> names;
	map<string,u64>::const_iterator it;
	for(it = routineMisses.begin(); it != routineMisses.end(); ++it) {
	    string n = it->first;
	    if (routineTotals[n] > 0)
		names.push_back(n);
	}

	sort(names.begin(), names.end(), namesorter(routineMisses));
	vector<string>::reverse_iterator rit;
	u64 missesShown = 0;
	for(rit = names.rbegin(); rit != names.rend(); ++rit) {
	    u64 total = routineTotals[*rit];
	    u64 misses = routineMisses[*rit];
	    u64 hits = total - misses;
	    fprintf(stderr, "  %30s : %8llu, %8llu misses (%6.2f %%), %8llu hits (%6.2f %%)\n",
		    (*rit).c_str(), total,
		    misses, 100.0 * misses/total,
		    hits, 100.0 * hits/total);

	    // stop after 99% misses shown
	    missesShown += misses;
	    if (100.0 * missesShown / threadMissCount > 99.0) break;
}

        for(int i=0; i<MAX_PREDCOUNTERS; i++) {
            if (d->predCounters[i] == 0) continue;
            fprintf(stderr, "   PredCounters %2d : %llu\n",
                    i, d->predCounters[i]);
        }
    }
    fprintf(stderr, "  Total predicted: %llu / %llu = %5.2f %%\n",
            totalHitCount, totalCount,
            100.0 * (double)totalHitCount / (double)totalCount);
}

/* ===================================================================== */
/* Main                                                                  */
/* ===================================================================== */

int main(int argc, char *argv[])
{
    if( PIN_Init(argc,argv) ) {
      fprintf(stderr, "PinPredict V0.1\n");
      return 0;
    }

    PIN_InitSymbols();
    predReg = PIN_ClaimToolRegister();
    storeFuncGen<MAX_STOREFUNC-1>();

    for(int i=0; i<MAX_THREADS; i++)
      tdata[i] = 0;

    TRACE_AddInstrumentFunction(Trace, 0);
    PIN_AddFiniFunction(Fini, 0);
    PIN_AddThreadStartFunction(ThreadStart, 0);

    // Never returns
    PIN_StartProgram();
    
    return 0;
}
