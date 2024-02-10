#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
/* #include <string.h> */
#include <sys/resource.h>
#include <sys/time.h>

#include "kseq.h"
#include "rlcsa.h"
#include "rle.h"
KSEQ_INIT(gzFile, gzread)

static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
    5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

double cputime() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
         1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + tp.tv_usec * 1e-6;
}

int main_2(int argc, char *argv[]) {
  char *fa_path = argv[1];
  rlcsa_t *rlc = rlc_init();

  gzFile fp = gzopen(fa_path, "rb");
  kseq_t *ks = kseq_init(fp);
  int64_t l;
  uint8_t *s;
  int i;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;
    printf("%s\n", s);

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
    for (i = 0; i < l; ++i)
      printf("%d", s[i]);
    printf("\n");

    rlc_insert(rlc, s, l);

    //   /* // Add reverse to buffer */
    //   /* for (i = 0; i < (l >> 1); ++i) { */
    //   /*   int tmp = s[l - 1 - i]; */
    //   /*   tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp; */
    //   /*   s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i]; */
    //   /*   s[i] = tmp; */
    //   /* } */
    //   /* if (l & 1) */
    //   /*   s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i]; */
    //   /* mr_insert1(mr, s); */
  }
  kseq_destroy(ks);
  gzclose(fp);

  rlc_destroy(rlc);

  return 0;
}

int main_1(int argc, char *argv[]) {
  int ret, i;
  uint8_t *block = calloc(0, sizeof(uint8_t));

  int64_t cnt[6];
  int64_t end_cnt[6] = {0, 0, 0, 0, 0, 0};
  int x = 0;
  int a = 2;
  int rl = 5;
  // insert the symbol $a after $x symbols in $block; marginal counts added to
  // $cnt; returns the size increase cnt are counts at insertion position
  // (changed during insertion) end_cnt are current counts at end (managed by
  // caller and increased after
  // insertion)
  ret = rle_insert(block, x, a, rl, cnt, end_cnt);
  printf("%d\n", ret);
  end_cnt[a] += rl;
  for (i = 0; i < 6; ++i)
    printf("%ld ", end_cnt[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    printf("%ld ", cnt[i]);
  printf("\n");

  x = rl - 1;
  a = 3;
  rl = 2;
  ret = rle_insert(block, x, a, rl, cnt, end_cnt);
  printf("%d\n", ret);
  end_cnt[a] += rl;
  for (i = 0; i < 6; ++i)
    printf("%ld ", end_cnt[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    printf("%ld ", cnt[i]);
  printf("\n");

  x = x + rl + 1;
  a = 4;
  rl = 4;
  ret = rle_insert(block, x, a, rl, cnt, end_cnt);
  printf("%d\n", ret);
  end_cnt[a] += rl;
  x += rl;
  for (i = 0; i < 6; ++i)
    printf("%ld ", end_cnt[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    printf("%ld ", cnt[i]);
  printf("\n\n");

  rle_print(block, 1);
  rle_print(block, 0);

  x = 0;
  int64_t cx[6] = {0, 0, 0, 0, 0, 0};
  rle_rank1a(block, x, cx, end_cnt);
  printf("%d: ", x);
  for (i = 0; i < 6; ++i)
    printf("%ld ", cx[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    cx[i] = 0;
  x = 1;
  rle_rank1a(block, x, cx, end_cnt);
  printf("%d: ", x);
  for (i = 0; i < 6; ++i)
    printf("%ld ", cx[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    cx[i] = 0;
  x = 3;
  rle_rank1a(block, x, cx, end_cnt);
  printf("%d: ", x);
  for (i = 0; i < 6; ++i)
    printf("%ld ", cx[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    cx[i] = 0;
  x = 4;
  rle_rank1a(block, x, cx, end_cnt);
  printf("%d: ", x);
  for (i = 0; i < 6; ++i)
    printf("%ld ", cx[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    cx[i] = 0;
  x = 9;
  rle_rank1a(block, x, cx, end_cnt);
  printf("%d: ", x);
  for (i = 0; i < 6; ++i)
    printf("%ld ", cx[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    cx[i] = 0;
  x = 10;
  rle_rank1a(block, x, cx, end_cnt);
  printf("%d: ", x);
  for (i = 0; i < 6; ++i)
    printf("%ld ", cx[i]);
  printf("\n");

  for (i = 0; i < 6; ++i)
    cx[i] = 0;
  x = 11;
  rle_rank1a(block, x, cx, end_cnt);
  printf("%d: ", x);
  for (i = 0; i < 6; ++i)
    printf("%ld ", cx[i]);
  printf("\n");

  return 0;
}

int main_0(int argc, char *argv[]) {
  double t_start;
  int i;
  uint64_t ret, l;
  t_start = realtime();
  srand(0);
  uint8_t *block = calloc((uint64_t)1 << 33, sizeof(uint8_t));
  int64_t cnt[6];
  int64_t end_cnt[6] = {0, 0, 0, 0, 0, 0};
  int p = 0;
  int a, rl;
  for (i = 0; i < 1000000; ++i) {
    a = rand() % 4 + 1;
    rl = rand() % 500 + 1;
    ret = rle_insert(block, p, a, rl, cnt, end_cnt);
    end_cnt[a] += rl;
    p += rl;
  }
  printf("%d\n", ret);
  for (i = 0; i < 6; ++i)
    printf("%ld ", end_cnt[i]);
  printf("\n");
  l = end_cnt[1] + end_cnt[2] + end_cnt[3] + end_cnt[4];
  fprintf(stderr, "%ld\n", l);
  // rle_print(block, 0);
  fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
          realtime() - t_start, cputime());

  // t_start = realtime();

  uint64_t cx[6] = {0, 0, 0, 0, 0, 0};
  rle_count(block, cx);
  for (i = 0; i < 6; ++i)
    printf("%ld ", cx[i]);
  printf("\n");
  l = cx[1] + cx[2] + cx[3] + cx[4];
  fprintf(stderr, "%ld\n", l);

  /* // for (i = 0; i < 6; ++i) */
  /* //   cx[i] = 0; */
  /* // for (int j = 0; j < 10000; ++j) */
  /* // { */
  /* //   p = rand() % l; */
  /* //   rle_rank1a(block, p, cx, end_cnt); */
  /* //   // printf("%d: ", p); */
  /* //   // for (i = 0; i < 6; ++i) */
  /* //   //   printf("%ld ", cx[i]); */

  /* //   // printf("\n"); */
  /* //   // for (i = 0; i < 6; ++i) */
  /* //   //   cx[i] = 0; */
  /* // } */

  free(block);
  /* // insert the symbol $a after $p symbols in $block; marginal counts added
   * to */
  /* // $cnt; returns the size increase cnt are counts at insertion position */
  /* // (changed during insertion) end_cnt are current counts at end (managed by
   */
  /* // caller and increased after insertion) */

  /* fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
   */
  /*         realtime() - t_start, cputime()); */

  return 0;
}

int main(int argc, char *argv[]) {
  main_0(argc - 1, argv + 1);
  return 0;
}
