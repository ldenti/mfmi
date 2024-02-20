#ifndef RLE6_H_
#define RLE6_H_

#include <stdint.h>
#include <stdio.h>

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x), 1)
#else
#define LIKELY(x) (x)
#endif

typedef struct {
  uint8_t *vector;
  uint8_t **blocks;
  uint32_t b;
  uint64_t l;
  uint32_t max_size;
  uint32_t sd;
  uint32_t *psamples;
  uint32_t *rsamples;
  int64_t *ecxs;
} rle_t;

/**
 * Initialize a rle vector
 *
 * @param size      maximum size for bit vectors
 * @param sd        sampling distance
 */
rle_t *rle_init(uint64_t size, int sd);

/**
 * Destroy the rle vector
 */
void rle_destroy(rle_t *rle);

/**
 * Insert a run in the vector
 *
 * @param rle       rle vector
 * @param x         position
 * @param a         character
 * @param rl        run length
 * @param cnt       counts at position x (updated by function)
 * @param end_cnt   counts at end of vector (managed by caller, must be
 * increased after insertion)
 */
int rle_insert(rle_t *rle, int64_t x, int a, int64_t rl, int64_t cnt[6],
               const int64_t end_cnt[6]);

/**
 * Build and sample ranks, freeze vector
 *
 * @param rle       rle vector
 */
void rle_sample(rle_t *rle);

/**
 * Get counts
 *
 * @param rle       rle vector
 * @param cnt       container for counts to be returned
 */
void rle_count(const rle_t *rle, int64_t cnt[6]);

// TODO: fix and improve this. For now it just works
void rle_rank2a(const uint8_t *vector, int64_t x, int64_t y, int64_t *cx,
                int64_t *cy, const int64_t ec[6]);

/**
 * Rank each symbol in [0..x)
 *
 * @param rle       rle vector
 * @param x         rank up to this position (excluded)
 * @param cx        container for ranks  to be returned
 */
void rle_rank1a(const rle_t *rle, int64_t x, int64_t *cx);

/**
 * Print the rle vector
 *
 * @param rle       rle vector
 * @param expand    expand (1) or not (0) the runs while printing
 */
void rle_print(const rle_t *rle, int expand);

/******************
 *** 43+3 codec ***
 ******************/

extern const uint8_t rle_auxtab[8];

#define RLE_MIN_SPACE 18
#define rle_nptr(vector) ((uint32_t *)(vector))

// decode one run (c,l) and move the pointer p
#define rle_dec1(p, c, l)                                                      \
  do {                                                                         \
    (c) = *(p) & 7;                                                            \
    if (LIKELY((*(p) & 0x80) == 0)) {                                          \
      (l) = *(p)++ >> 3;                                                       \
    } else if (LIKELY(*(p) >> 5 == 6)) {                                       \
      (l) = (*(p) & 0x18L) << 3L | ((p)[1] & 0x3fL);                           \
      (p) += 2;                                                                \
    } else {                                                                   \
      int n = ((*(p) & 0x10) >> 2) + 4;                                        \
      (l) = *(p)++ >> 3 & 1;                                                   \
      while (--n)                                                              \
        (l) = ((l) << 6) | (*(p)++ & 0x3fL);                                   \
    }                                                                          \
  } while (0)

static inline int rle_enc1(uint8_t *p, int c, int64_t l) {
  if (l < 1LL << 4) { // 16 -> 1 byte
    *p = l << 3 | c;
    return 1;
  } else if (l < 1LL << 8) { // 256 -> 2 bytes
    *p = 0xC0 | l >> 6 << 3 | c;
    p[1] = 0x80 | (l & 0x3f);
    return 2;
  } else if (l < 1LL << 19) { // 524228 -> 4 bytes
    *p = 0xE0 | l >> 18 << 3 | c;
    p[1] = 0x80 | (l >> 12 & 0x3f);
    p[2] = 0x80 | (l >> 6 & 0x3f);
    p[3] = 0x80 | (l & 0x3f);
    return 4;
  } else { // 8 bytes
    int i, shift = 36;
    *p = 0xF0 | l >> 42 << 3 | c;
    for (i = 1; i < 8; ++i, shift -= 6)
      p[i] = 0x80 | (l >> shift & 0x3f);
    return 8;
  }
}

#endif
