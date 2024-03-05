#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"

KSEQ_INIT(gzFile, gzread)

void pp(rlcsa_t *rlc, const uint8_t *seq, const int64_t l, const char *qname) {
  int begin = l - 1;
  uint8_t first = 1;
  uint8_t c;
  bisa_t sai;
  while (begin >= 0) {
    c = seq[begin];
    rlc_init_biinterval(rlc, c, sai);
    fprintf(stderr, "Starting from %c (%d): %d,%d,%d\n", "$ACGTN"[c], begin,
            sai.x[0], sai.x[1], sai.x[2]);

    // Backward search. Stop at first mismatch.
    int bmatches = 0;
    while (sai.x[2] > 0 && begin > 0) {
      --begin;
      c = seq[begin];
      ++bmatches;
      bisa_t osai[6];
      rlc_bilf(rlc, &sai, osai, 1);
      sai = osai[seq[begin]];
      fprintf(stderr, "Backward extending with %c (%d): %d,%d,%d\n",
              "$ACGTN"[c], begin, sai.x[0], sai.x[1], sai.x[2]);
      // fprintf(stderr, "%d,%d,%d\n", sai.x, sai.rx, sai.l);
    }
    // prefix is a match
    if (begin == 0 && sai.x[2] > 0)
      break;

    // Forward search
    int end = begin;
    int fmatches = 0;
    c = seq[end];
    rlc_init_biinterval(rlc, c, sai);
    fprintf(stderr, "Starting from %c (%d): %d,%d,%d\n", "$ACGTN"[c], end,
            sai.x[0], sai.x[1], sai.x[2]);
    while (sai.x[2] > 0) {
      ++end;
      c = seq[end];
      ++fmatches;
      bisa_t osai[6];
      rlc_bilf(rlc, &sai, osai, 0);
      sai = osai[fm6_comp(seq[end])];
      fprintf(stderr, "Forward extending with %c (%d): %d,%d,%d\n", "$ACGTN"[c],
              end, sai.x[0], sai.x[1], sai.x[2]);
    }

    fprintf(stderr, "%d -> %d\n", begin, end);
    fprintf(stderr, "Matches: <%d | >%d\n", bmatches, fmatches);

    // add solution
    printf("%s\t%d\t%d\n", first ? qname : "*", begin, end - begin + 1);
    /* for (int i = begin; i < end; ++i) */
    /*   printf("%c", "$ACGTN"[seq[i]]); */
    /* printf("\n"); */
    first = 0;

    // TODO: improve this

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
  char *fq_fn = argv[2];    // **perfect** reads

  rlcsa_t *rlc = rlc_restore(index_fn);

  gzFile fp = gzopen(fq_fn, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  uint8_t *s;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

    pp(rlc, s, l, ks->name.s);
  }
  kseq_destroy(ks);
  gzclose(fp);
  rlc_destroy(rlc);

  fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
          realtime() - t_start, cputime());

  return 0;
}
