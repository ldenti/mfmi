#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.hpp"

KSEQ_INIT(gzFile, gzread)

void pp(rlcsa_t *rlc, const uint8_t *seq, const int64_t l, const char *qname) {
  int begin = l - 1;
  uint8_t first = 1;
  uint8_t c;
  qint_t sai;
  qint_t osai[6];
  while (begin >= 0) {
    c = seq[begin];
    rlc_init_qinterval(rlc, c, sai);
    // fprintf(stderr, "Starting from %c (%d): %ld,%ld,%ld\n", "$ACGTN"[c],
    // begin, sai.x[0], sai.x[1], sai.x[2]);

    // Backward search. Stop at first mismatch.
    while (sai.x[2] > 0 && begin > 0) {
      --begin;
      c = seq[begin];
      rlc_extend(rlc, &sai, osai, 1);
      sai = osai[seq[begin]];
      // fprintf(stderr, "Backward extending with %c (%d): %ld,%ld,%ld\n",
      //         "$ACGTN"[c], begin, sai.x[0], sai.x[1], sai.x[2]);
    }
    // prefix is a match
    if (begin == 0 && sai.x[2] > 0)
      break;

    // Forward search
    int end = begin;
    c = seq[end];
    rlc_init_qinterval(rlc, c, sai);
    // fprintf(stderr, "Starting from %c (%d): %ld,%ld,%ld\n", "$ACGTN"[c],
    // end,
    //         sai.x[0], sai.x[1], sai.x[2]);
    while (sai.x[2] > 0) {
      ++end;
      c = seq[end];
      rlc_extend(rlc, &sai, osai, 0);
      sai = osai[fm6_comp(seq[end])];
      // fprintf(stderr, "Forward extending with %c (%d): %ld,%ld,%ld\n",
      //         "$ACGTN"[c], end, sai.x[0], sai.x[1], sai.x[2]);
    }

    // add solution
    printf("%s\t%d\t%d\n", first ? qname : "*", begin, end - begin + 1);
    // for (int i = begin; i < end; ++i)
    //   printf("%c", "$ACGTN"[seq[i]]);
    // printf("\n");
    first = 0;

    // prepare for next round
    if (begin == 0)
      break;
    begin = end - 1;
    // TODO: add relaxed and overlap (?)
  }
}

int main_pingpong(int argc, char *argv[]) {
  (void)argc; // suppress unused parameter warning
  double t_start;
  t_start = realtime();

  char *index_fn = argv[1]; // index
  char *fq_fn = argv[2];    // queries

  rlcsa_t *rlc = rlc_restore(index_fn);

  gzFile fp = gzopen(fq_fn, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  uint8_t *s;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = fm6(s[i]);

    pp(rlc, s, l, ks->name.s);
  }
  kseq_destroy(ks);
  gzclose(fp);
  rlc_destroy(rlc);

  fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
          realtime() - t_start, cputime());

  return 0;
}
