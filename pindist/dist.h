/* Calculation of a Stack Reuse Distance Histogram
 *
 * (C) 2015, Josef Weidendorfer / LRR-TUM
 * GPLv2+ (see COPYING)
 */

typedef void* Addr;

// initialize / clear used structs
void RD_init(int min1);

// add distance buckets, starting from smallest (>0)
// only specification of minimal distance required
// use as last bucket min=0 to get "infinite distances" there
void RD_addBucket(unsigned int min);

// run stack simulation
// To use a specific block size, ensure that <a> is aligned
void RD_accessBlock(Addr a);

// do an internal consistency check
void RD_checkConsistency();

// get statistics
void RD_stat(unsigned long & stack_size, unsigned long & accessCount);

// get resulting histogram
// Repeatly call RD_get_hist, start with bucket 0.
// Returns next bucket or 0 if this was last
int RD_get_hist(unsigned int bucket,
		unsigned int & min, unsigned long & accessCount);

// print nice ASCII histogram to <out>
//  <pStr> is prefix for every line, distances scaled by <blockSize>
void RD_printHistogram(FILE* out, const char* pStr, int blockSize);
