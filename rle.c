#include "rle.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const uint8_t rle_auxtab[8] = {0x01, 0x11, 0x21, 0x31, 0x03, 0x13, 0x07, 0x17};

rle_t *rle_init(uint64_t size, int sd) {
  rle_t *rle = (rle_t *)calloc(1, sizeof(rle_t));
  rle->max_size = size;
  rle->sd = sd;
  rle->b = 0;
  rle->l = 0;
  rle->vector = (uint8_t *)calloc(size + 4, sizeof(uint8_t));
  rle->blocks = (uint8_t **)calloc(size / sd + 1, sizeof(uint8_t *));
  rle->psamples = (uint32_t *)calloc(size / sd + 1, sizeof(uint32_t));
  rle->rsamples = (uint32_t *)calloc(size / sd + 1, sizeof(uint32_t));
  rle->ecxs = (int64_t *)calloc(6 * (size / sd + 1), sizeof(int64_t));
  return rle;
}

void rle_destroy(rle_t *rle) {
  free(rle->vector);
  free(rle->blocks);
  free(rle->psamples);
  free(rle->rsamples);
  free(rle->ecxs);
  free(rle);
}

// insert symbol $a after $x symbols in $str; marginal counts added to $cnt;
// returns the size increase
int rle_insert_cached(uint8_t *vector, int64_t x, int a, int64_t rl,
                      int64_t cnt[6], const int64_t ec[6], int *beg,
                      int64_t bc[6]) {
  uint32_t *nptr = (uint32_t *)vector;
  int diff;

  vector += 4; // skip the first 4 counting bytes
  if (*nptr == 0) {
    memset(cnt, 0, 48);
    diff = rle_enc1(vector, a, rl);
  } else {
    uint8_t *p, *end = vector + *nptr, *q;
    int64_t pre, z, l = 0, tot, beg_l;
    int c = -1, n_bytes = 0, n_bytes2, t = 0;
    uint8_t tmp[24];
    beg_l = bc[0] + bc[1] + bc[2] + bc[3] + bc[4] + bc[5];
    tot = ec[0] + ec[1] + ec[2] + ec[3] + ec[4] + ec[5];
    if (x < beg_l) {
      beg_l = 0, *beg = 0;
      memset(bc, 0, 48);
    }
    if (x == beg_l) {
      p = q = vector + (*beg);
      z = beg_l;
      memcpy(cnt, bc, 48);
    } else if (x - beg_l <=
               ((tot - beg_l) >> 1) + ((tot - beg_l) >> 3)) { // forward
      z = beg_l;
      p = vector + (*beg);
      memcpy(cnt, bc, 48);
      while (z < x) {
        rle_dec1(p, c, l);
        z += l;
        cnt[c] += l;
      }
      for (q = p - 1; *q >> 6 == 2; --q)
        ;
    } else { // backward
      memcpy(cnt, ec, 48);
      z = tot;
      p = end;
      while (z >= x) {
        --p;
        if (*p >> 6 != 2) {
          l |= *p >> 7 ? (int64_t)rle_auxtab[*p >> 3 & 7] >> 4 << t : *p >> 3;
          z -= l;
          cnt[*p & 7] -= l;
          l = 0;
          t = 0;
        } else {
          l |= (*p & 0x3fL) << t;
          t += 6;
        }
      }
      q = p;
      rle_dec1(p, c, l);
      z += l;
      cnt[c] += l;
    }
    *beg = q - vector;
    memcpy(bc, cnt, 48);
    bc[c] -= l;
    n_bytes = p - q;
    if (x == z && a != c && p < end) { // then try the next run
      int tc;
      int64_t tl;
      q = p;
      rle_dec1(q, tc, tl);
      if (a == tc)
        c = tc, n_bytes = q - p, l = tl, z += l, p = q, cnt[tc] += tl;
    }
    if (z != x)
      cnt[c] -= z - x;
    pre = x - (z - l);
    p -= n_bytes;
    if (a == c) { // insert to the same run
      n_bytes2 = rle_enc1(tmp, c, l + rl);
    } else if (x == z) { // at the end; append to the existing run
      p += n_bytes;
      n_bytes = 0;
      n_bytes2 = rle_enc1(tmp, a, rl);
    } else { // break the current run
      n_bytes2 = rle_enc1(tmp, c, pre);
      n_bytes2 += rle_enc1(tmp + n_bytes2, a, rl);
      n_bytes2 += rle_enc1(tmp + n_bytes2, c, l - pre);
    }
    if (n_bytes != n_bytes2 && end != p + n_bytes) // size changed
      memmove(p + n_bytes2, p + n_bytes, end - p - n_bytes);
    memcpy(p, tmp, n_bytes2);
    diff = n_bytes2 - n_bytes;
  }
  return (*nptr += diff);
}

int rle_insert(rle_t *rle, int64_t x, int a, int64_t rl, int64_t cnt[6],
               const int64_t ec[6]) {
  int beg = 0;
  int64_t bc[6];
  memset(bc, 0, 48);
  return rle_insert_cached(rle->vector, x, a, rl, cnt, ec, &beg, bc);
}

void rle_update_block(rle_t *rle, uint8_t c, uint32_t l) {
  // printf("Updating block %d, character %d (p=%d): %d ", rle->b, c, rle->b * 6
  // + c, rle->ecxs[(rle->b * 6 + c)]);
  rle->ecxs[(rle->b * 6 + c)] += l;
  // printf("-> %d\n", rle->ecxs[(rle->b * 6 + c)]);
}

void rle_store_block(rle_t *rle, uint8_t *q, uint32_t p, uint32_t r) {
  rle->blocks[rle->b] = q;
  rle->psamples[rle->b] = p;
  rle->rsamples[rle->b] = r;
  ++rle->b;
}

void rle_sample(rle_t *rle) {
  const uint32_t *n = (const uint32_t *)rle->vector;
  uint8_t *q = rle->vector + 4, *end = rle->vector + 4 + *n;
  uint8_t *pq = q;
  int i = 0;          // each run can take more bytes, this counts the bytes
  int rank = 0;       // sample for rank
  int currp = 0;      // sample for positions
  int currb_rank = 0; // rank in current block
  uint8_t c;
  uint32_t l;
  while (q < end) {
    if (i >= rle->sd) {
      rle_store_block(rle, pq, currp, rank);
      i = 0;
      currp = rle->l; // += q - pq;
      pq = q;
      rank += currb_rank;
      currb_rank = 0;
    }
    rle_dec1(q, c, l);
    rle->l += l;
    if (c)
      currb_rank += l;
    rle_update_block(rle, c, l);
    ++i;
  }
  rle_store_block(rle, pq, currp, rank);
  /* printf("%d BLOCKS\n", rle->b); */
  /* for (int b = 0; b < rle->b; ++b) */
  /*   rle_print_block(rle, b, 1); */
}

void rle_count(const rle_t *rle, int64_t cnt[6]) {
  const uint8_t *q = rle->vector + 4, *end = q + *(uint32_t *)rle->vector;
  while (q < end) {
    int c;
    int64_t l;
    rle_dec1(q, c, l);
    cnt[c] += l;
  }
}

void rle_print(const rle_t *rle, int expand) {
  const uint32_t *p = (const uint32_t *)rle->vector;
  const uint8_t *q = rle->vector + 4, *end = rle->vector + 4 + *p;
  while (q < end) {
    int c;
    int64_t l, x;
    rle_dec1(q, c, l);
    if (expand)
      for (x = 0; x < l; ++x)
        putchar("01XXXX"[c]);
    else
      printf("%c(%ld) ", "01XXXX"[c], (long)l);
  }
  putchar('\n');
}

void rle_rank1a(const rle_t *rle, int64_t x, int64_t *cx) {
  // rle_print(rle, 1);
  // int _b = x / rle->sd; // block index
  int b = -1;
  for (int i = 0; i < rle->b; ++i) {
    // printf("%d %d %d\n", rle->psamples[i], x, rle->psamples[i] < x);
    b += rle->psamples[i] < x;
  }
  // printf("b: %d\n", b);
  /* for (int i = 0; i < rle->b; ++i) */
  /*   printf("ranks %d: %d\n", i, rle->rsamples[i]); */
  /* for (int i = 0; i < rle->b; ++i) */
  /*   printf("poss %d: %d\n", i, rle->psamples[i]); */

  /* printf("EXC[0]: %d\n", rle->ecxs[(b * 6)]); */
  /* printf("EXC[1]: %d\n", rle->ecxs[(b * 6 + 1)]); */
  /* if (x >= rle->l) { */
  /*   cx[0] = rle->ecxs[b * 6 + 0]; */
  /*   cx[1] = rle->ecxs[b * 6 + 1]; */
  /* } else { */
  rle_rank2a(rle->blocks[b], x - rle->psamples[b], -1, cx, 0,
             (rle->ecxs + (b * 6)));
  /* printf("Global p: %d\n", x); */
  /* printf("  - %d\n", rle->psamples[b]); */
  /* printf("Local p: %d\n", x - rle->psamples[b]); */

  /* printf("Local rank: %d\n", cx[1]); */
  /* printf("  + %d\n", rle->rsamples[b]); */
  cx[1] += rle->rsamples[b];
  /* } */
}

void rle_rank2a(const uint8_t *vector, int64_t x, int64_t y, int64_t *cx,
                int64_t *cy, const int64_t ec[6]) {
  /* printf("'' EXC[0]: %d\n", ec[0]); */
  /* printf("'' EXC[1]: %d\n", ec[1]); */
  int a;
  int64_t tot, cnt[6];
  const uint8_t *p;

  y = y >= x ? y : x;
  tot = ec[0] + ec[1] + ec[2] + ec[3] + ec[4] + ec[5];
  if (tot == 0)
    return;
  // printf("RANK2 CHECK: %d <= (%d - %d) + %d\n", x, tot, y, tot >> 3);
  if (1) { // x <= (tot - y) + (tot >> 3)) {
    int c = 0;
    int64_t l, z = 0;
    memset(cnt, 0, 48);
    p = vector; //  + 4;
    while (z < x) {
      rle_dec1(p, c, l);
      z += l;
      cnt[c] += l;
    }
    for (a = 0; a != 6; ++a)
      cx[a] += cnt[a];
    cx[c] -= z - x;
    if (cy) {
      while (z < y) {
        rle_dec1(p, c, l);
        z += l;
        cnt[c] += l;
      }
      for (a = 0; a != 6; ++a)
        cy[a] += cnt[a];
      cy[c] -= z - y;
    }
  } else {
#define move_backward(_x)                                                      \
  while (z >= (_x)) {                                                          \
    --p;                                                                       \
    if (*p >> 6 != 2) {                                                        \
      l |= *p >> 7 ? (int64_t)rle_auxtab[*p >> 3 & 7] >> 4 << t : *p >> 3;     \
      z -= l;                                                                  \
      cnt[*p & 7] -= l;                                                        \
      l = 0;                                                                   \
      t = 0;                                                                   \
    } else {                                                                   \
      l |= (*p & 0x3fL) << t;                                                  \
      t += 6;                                                                  \
    }                                                                          \
  }

    int t = 0;
    int64_t l = 0, z = tot;
    memcpy(cnt, ec, 48);
    p = vector; // + *(const uint32_t *)vector; // +4
    if (cy) {
      move_backward(y) for (a = 0; a != 6; ++a) cy[a] += cnt[a];
      cy[*p & 7] += y - z;
    }
    move_backward(x) for (a = 0; a != 6; ++a) cx[a] += cnt[a];
    cx[*p & 7] += x - z;

#undef move_backward
  }
}
