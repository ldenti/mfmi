#ifndef RLCSA_H_
#define RLCSA_H_

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include "ksort.h"
#include "kvec.h"
#include "rle.h"
#include "rope.h"
#include "utils.h"

// void radix_sort(int64_t *array, int64_t offset, int64_t end, int shift);

typedef struct {
  int64_t a;
  int64_t b;
} pair_t;

typedef struct {
  int64_t x;
  int64_t rx;
  int64_t l;
} bisa_t;

#define pair_lt(x, y) ((x).b < (y).b)

// FIXME: do we want these?
typedef pair_t sa_t;
typedef pair_t skew_pair;
typedef pair_t ss_range;

// structs to allow kvec to be passed around
// (see https://github.com/attractivechaos/klib/issues/144)
typedef struct ss_ranges {
  kvec_t(ss_range);
} ss_ranges;

typedef struct uint_kv {
  kvec_t(int64_t);
} int_kv;

typedef struct {
  int64_t l; // length of bitvectors in bits (same as length of indexed text)
  int64_t *cnts;   // counts for each symbol
  int64_t *C;      // C array
  rope_t *bits[6]; // the bit vectors
} rlcsa_t;         // the rlcsa index

/**
 * Initialize a rlcsa
 */
rlcsa_t *rlc_init();

/**
 * Destroy the rlcsa index
 */
void rlc_destroy(rlcsa_t *rlc);

/**
 * Write rlcsa index to file
 */
int rlc_dump(rlcsa_t *rlc, const char *fn);

/**
 * Restore rlcsa index from file. Index must be freed by caller.
 */
rlcsa_t *rlc_restore(const char *fn);

/**
 * Insert multiple (0-separated) strings into the index
 *
 * @param rlc       the index
 * @param seq       the input string(s)
 * @param n         input string length (last 0 not included)
 */
void rlc_insert(rlcsa_t *rlc, const uint8_t *seq, int64_t n);

/**
 * Return q-interval for a character
 *
 * @param rlc       the index
 * @param c         the character
 */
sa_t rlc_init_interval(rlcsa_t *rlc, uint8_t c);

/**
 * LF-mapping (backward extension)
 *
 * @param rlc       the index
 * @param range     a q-interval
 * @param c         a character
 */
sa_t rlc_lf(rlcsa_t *rlc, sa_t range, uint8_t c);

/**
 * Return bi-interval for a character (FMD)
 *
 * @param rlc       the index
 * @param c         the character
 */
bisa_t rlc_init_biinterval(rlcsa_t *rlc, uint8_t c);

/**
 * LF-mapping (backward extension) for bi-intervals (FMD)
 *
 * @param rlc       the index
 * @param range     a q-interval
 * @param c         the character
 * @param backward  backward/forward extension
 */
bisa_t rlc_bilf(rlcsa_t *rlc, bisa_t range, uint8_t c, uint8_t backward);

/**
 * Merge rlc2 into rlc1. rlc2 is freed
 *
 * @param rlc1      first index
 * @param rlc2      second index
 * @param seq       the text of rlc2
 */
void rlc_merge(rlcsa_t *rlc, rlcsa_t *increment, const uint8_t *seq);

#endif
