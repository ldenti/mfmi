#ifndef RLCSA_H_
#define RLCSA_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "kvec.h"
#include "rle.h"

typedef struct {
  uint32_t a;
  uint32_t b;
} sa_t;

// CHECKME: do we want these?
typedef sa_t skew_pair;
typedef sa_t ss_range;

// struct to allow kvec to be passed around
// (see https://github.com/attractivechaos/klib/issues/144)
typedef struct ss_ranges {
  kvec_t(ss_range);
} ss_ranges;

typedef struct {
  uint64_t *cnts;
  uint64_t *C;
  rle_t *bits[6];
  int64_t *end_cnts[6];
} rlcsa_t;

/**
 * Initialize a rlcsa
 *
 * @param size      maximum size for bit vectors
 * @param sd        sampling distance
 */
rlcsa_t *rlc_init(uint32_t size, uint32_t sd);

/**
 * Insert one string into the index
 *
 * @param rlc       the index
 * @param seq       the input string (no $-terminated)
 * @param n         input string length
 */
int rlc_insert(rlcsa_t *rlc, const uint8_t *seq, int n);

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
sa_t LF(rlcsa_t *rlc, sa_t range, uint8_t c);

/**
 * Destroy the rlcsa index
 */
void rlc_destroy(rlcsa_t *rlc);

#endif
