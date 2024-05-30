#ifndef RLCSA_HPP_
#define RLCSA_HPP_

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <execution>
#include <omp.h>
#include <vector>

#include "bits/nibblevector.h"
#include "rld0.h"
#include "utils.h"

// #include <assert.h>

// #include <stdint.h>
// #include <stdlib.h>

// FIXME: do we want these?
typedef std::pair<uint32_t, uint32_t> sa_t;
typedef std::pair<uint32_t, uint32_t> skew_pair;
typedef std::pair<uint32_t, uint32_t> ss_range;

// typedef NibbleVector  PsiVector;
// typedef RLEVector     PsiVector;

typedef struct { // the rlcsa index
  int64_t l;     // length of indexed text
  int64_t *cnts; // marginal counts for each symbol (same as rope->c)
  int64_t *C;    // C array
  CSA::NibbleVector *bits[6]; // the bit vectors representing the BWT
} rlcsa_t;

/**
 * Initialize a rlcsa
 */
rlcsa_t *rlc_init();

/**
 * Destroy the rlcsa index
 */
void rlc_destroy(rlcsa_t *rlc);

/**
 * Dump rlc to file fp
 */
void rlc_dump(rlcsa_t *rlc, const char *fp);

/**
 * Print bwt to stdout
 */
void rlc_print_bwt(rlcsa_t *rlc);

/**
 * Insert multiple (0-separated) strings into the index. If index is empty,
 * build it. Otherwise, build a new index and merge.
 *
 * @param rlc       index to update
 * @param seq       the input string(s)
 * @param n         input string length (last 0 not included)
 * @param nt        number of threads
 */
void rlc_insert(rlcsa_t *rlc, const uint8_t *seq, uint32_t n, int nt);

/**
 * Build index from sequence.
 *
 * @params rlc       the index
 * @params sequence  the sequence
 * @params n         sequence length
 * @params nt        number of threads
 */
void rlc_build(rlcsa_t *rlc, const uint8_t *sequence, uint32_t n, int nt);

/**
 * Merge rlc2 into rlc1. rlc2 is freed
 *
 * @param rlc1      first index
 * @param rlc2      second index
 * @param seq       the text of rlc2
 * @param nt        number of threads
 */
void rlc_merge(rlcsa_t *rlc, rlcsa_t *increment, const uint8_t *seq, int nt);

#endif
