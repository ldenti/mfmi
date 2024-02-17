#ifndef RLCSA_H_
#define RLCSA_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// #include "ksort.h"
#include "kvec.h"
#include "rle.h"

typedef struct {
  uint32_t a;
  uint32_t b;
} sa_t;

typedef sa_t skew_pair;
typedef sa_t ss_range;

typedef struct ss_ranges {
  kvec_t(ss_range);
} ss_ranges;

typedef struct {
  uint64_t *cnts;
  uint64_t *C;
  rle_t *bits[6];
  int64_t *end_cnts[6];
} rlcsa_t;

rlcsa_t *rlc_init();

int rlc_insert(rlcsa_t *rlc, const uint8_t *seq, int n);

sa_t rlc_init_interval(rlcsa_t *rlc, uint8_t c);
sa_t LF(rlcsa_t *rlc, sa_t range, uint8_t c);

sa_t *simpleSuffixSort(const uint8_t *sequence, uint32_t n, uint8_t threads);

void rlc_destroy(rlcsa_t *rlc);

#endif
