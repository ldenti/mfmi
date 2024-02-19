#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "rle.h"

#include "kvec.h"

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

int main_testkvec(int argc, char *argv[]) {

  kvec_t(int) v;
  kv_push(int, v, 0);
  kv_push(int, v, 1);
  kv_push(int, v, 2);

  for (uint i = 0; i < kv_size(v); ++i) {
    printf("%d ", kv_A(v, i));
  }
  printf("\n");
  printf("%d %d\n", kv_size(v), v.n);
  v.n = 0;
  for (uint i = 0; i < kv_size(v); ++i) {
    printf("%d ", kv_A(v, i));
  }
  printf("\n");
  printf("%d %d\n", kv_size(v), v.n);

  kv_push(int, v, 4);
  kv_push(int, v, 5);

  for (uint i = 0; i < kv_size(v); ++i) {
    printf("%d ", kv_A(v, i));
  }
  printf("\n");
  printf("%d %d\n", kv_size(v), v.n);
  kv_destroy(v);

  return 0;
}

int main(int argc, char *argv[]) {
  double t_start;
  t_start = realtime();

  char *fa_path = argv[1]; // single reference
  char *fq_path = argv[2]; // reads on + strand

  rlcsa_t *rlc = rlc_init((uint32_t)1 << 31, 512);

  gzFile fp = gzopen(fa_path, "rb");
  kseq_t *ks = kseq_init(fp);
  int64_t l;
  uint8_t *s;
  int i;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
    rlc_insert(rlc, (const uint8_t *)s, l);

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

  fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
          realtime() - t_start, cputime());
  t_start = realtime();

  for (i = 1; i < 5; ++i) {
    // rle_print(rlc->bits[1], 0);
    // rle_print(rlc->bits[1], 1);
    rle_freeze(rlc->bits[i]);
    for (int b = 0; b < rlc->bits[i]->b; ++b) {
      // printf("--- %d\n", b);
      /* printf("-\n", b); */
      /* for (int j = 0; j < 6; ++j) */
      /*   printf("%d ", rlc->bits[i]->ecxs[b * 6 + j]); */
      /* printf("\n"); */
    }
  }

  sa_t interval;
  for (int c = 0; c < 6; ++c) {
    interval = rlc_init_interval(rlc, c);
    // printf("%d: [%d, %d]\n", c, interval.a, interval.b);
  }

  int hit;
  fp = gzopen(fq_path, "rb");
  ks = kseq_init(fp);
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;
    // printf("%s\n", ks->name.s);
    // printf("%s\n", ks->seq.s);

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
    i = l - 1;
    interval = rlc_init_interval(rlc, s[i]);
    // fprintf(stderr, "(1) %d: %d %d\n", s[i], interval.a, interval.b);
    --i;
    hit = 1;

    for (; i >= 0; --i) {
      interval = LF(rlc, interval, s[i]);
      // fprintf(stderr, "(2) %d: %d %d\n", s[i], interval.a, interval.b);
      if (interval.b < interval.a) {
        hit = 0;
        break;
      }
    }

    printf("%d\n", hit);
  }

  rlc_destroy(rlc);

  fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
          realtime() - t_start, cputime());

  return 0;
}

int main_debug(int argc, char *argv[]) {
  char *fa_path = "example/tiny.fa";

  rlcsa_t *rlc = rlc_init(1024, 2);

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
    printf("0\n");
    rlc_insert(rlc, (const uint8_t *)s, l);
  }
  kseq_destroy(ks);
  gzclose(fp);

  for (i = 1; i < 5; ++i) {
    printf("----- %d\n", i);
    rle_print(rlc->bits[i], 0);
    rle_freeze(rlc->bits[i]);
    for (int b = 0; b < rlc->bits[i]->b; ++b) {
      printf("--- %d\n", b);
      for (int j = 0; j < 6; ++j)
        printf("%d ", rlc->bits[i]->ecxs[b * 6 + j]);
      printf("\n");
    }
    printf("---\n");
  }

  sa_t interval;
  for (int c = 0; c < 6; ++c) {
    interval = rlc_init_interval(rlc, c);
    printf("%d: [%d, %d]\n", c, interval.a, interval.b);
  }

  char *pattern = "GATAT";
  printf("%s\n", pattern);
  l = strlen(pattern);
  i = l - 1;
  interval = rlc_init_interval(rlc, seq_nt6_table[pattern[i]]);
  printf("%d: [%d, %d]\n", seq_nt6_table[pattern[i]], interval.a, interval.b);
  --i;
  for (; i >= 0; --i) {
    printf("\n\n# i = %d, %c (%d)\n", i, pattern[i], seq_nt6_table[pattern[i]]);
    interval = LF(rlc, interval, seq_nt6_table[pattern[i]]);
    if (interval.b < interval.a) {
      printf("XXX %d: [%d, %d]\n", seq_nt6_table[pattern[i]], interval.a,
             interval.b);
      break;
    }
    printf("%d: [%d, %d]\n", seq_nt6_table[pattern[i]], interval.a, interval.b);
  }
  rlc_destroy(rlc);

  return 0;
}
