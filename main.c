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

static inline uint kputsn(const char *p, uint l, kstring_t *s) {
  if (s->l + l + 1 >= s->m) {
    char *tmp;
    s->m = s->l + l + 2;
    kroundup32(s->m);
    if ((tmp = (char *)realloc(s->s, s->m)))
      s->s = tmp;
    else
      return EOF;
  }
  memcpy(s->s + s->l, p, l);
  s->l += l;
  s->s[s->l] = 0;
  return l;
}

double cputime() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return (double)r.ru_utime.tv_sec + (double)r.ru_stime.tv_sec +
         1e-6 * (double)(r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return (double)tp.tv_sec + (double)tp.tv_usec * 1e-6;
}

int main(int argc, char *argv[]) {
  double t_start;
  t_start = realtime();

  char *fa_path = argv[1]; // reference
  char *fq_path = argv[2]; // reads on + strand

  rlcsa_t *rlc = rlc_init((uint64_t)1 << 32, 512);

  gzFile fp = gzopen(fa_path, "rb");
  kseq_t *ks = kseq_init(fp);
  int l;
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  int i;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

    // Add forward to buffer
    kputsn((char *)ks->seq.s, ks->seq.l + 1, &buf);

    // Add reverse to buffer
    for (i = 0; i < (l >> 1); ++i) {
      int tmp = s[l - 1 - i];
      tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
      s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
      s[i] = tmp;
    }
    if (l & 1)
      s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
    kputsn((char *)s, ks->seq.l + 1, &buf);
  }

//  printf("%d\n", buf.l);
//  for (i = 0; i < buf.l + 1; ++i)
//    printf("%d", buf.s[i]);
//  printf("\n");

  rlc_insert(rlc, (const uint8_t *)buf.s, buf.l);

  free(buf.s);
  kseq_destroy(ks);
  gzclose(fp);

  fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
          realtime() - t_start, cputime());
  t_start = realtime();

  for (i = 1; i < 5; ++i) {
    // rle_print(rlc->bits[1], 0);
    // rle_print(rlc->bits[1], 1);
    rle_sample(rlc->bits[i]);
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
