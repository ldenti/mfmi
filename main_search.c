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

  qint_t interval;
  for (int c = 0; c < 6; ++c) {
    rlc_init_qinterval(rlc, c, interval);
    printf("%d: [%ld, %ld]\n", c, interval.x[0], interval.x[0] + interval.x[2]);
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
    rlc_init_qinterval(rlc, s[i], interval);
    --i;
    for (; i >= 0; --i) {
      qint_t osai[6];
      rlc_extend(rlc, &interval, osai, 1);
      interval = osai[s[i]];
      if (interval.x[2] <= 0) {
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
