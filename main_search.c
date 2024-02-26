#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"

KSEQ_INIT(gzFile, gzread)

int main_search(int argc, char *argv[]) {
  (void)argc; // suppress unused parameter warning
  double t_start;
  t_start = realtime();

  char *index_fn = argv[1]; // index
  char *fq_fn = argv[2];    // **perfect** reads

  rlcsa_t *rlc = rlc_restore(index_fn);

  sa_t interval;
  for (int c = 0; c < 6; ++c) {
    interval = rlc_init_interval(rlc, c);
    printf("%d: [%ld, %ld]\n", c, interval.a, interval.b);
  }

  int errors = 0;
  gzFile fp = gzopen(fq_fn, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  uint8_t *s;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
    i = l - 1;
    interval = rlc_init_interval(rlc, s[i]);
    --i;
    for (; i >= 0; --i) {
      interval = rlc_lf(rlc, interval, s[i]);
      if (interval.b < interval.a) {
        errors += 1;
        break;
      }
    }
  }
  printf("Errors: %d\n", errors);
  kseq_destroy(ks);
  gzclose(fp);
  rlc_destroy(rlc);

  fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
          realtime() - t_start, cputime());

  return 0;
}
