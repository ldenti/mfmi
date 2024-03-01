#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "rld0.h"

#include "utils.h"

KSEQ_INIT(gzFile, gzread)

void pp(rld_t *index, const uint8_t *seq, const int64_t l, const char *qname) {
  int first = 1;
  rldintv_t sai;
  int begin = l - 1;
  while (begin >= 0) {
    // Backward search. Stop at first mismatch.
    fm6_set_intv(index, seq[begin], sai);
    while (sai.x[2] != 0 && begin > 0) {
      begin--;
      rldintv_t
          osai[6]; // output SA intervals (one for each symbol between 0 and 5)
      rld_extend(index, &sai, osai, 1);
      sai = osai[seq[begin]];
    }
    // last checked char (i.e., first of the query) was a match
    if (begin == 0 && sai.x[2] != 0)
      break;

    //  Forward search:
    int end = begin;
    fm6_set_intv(index, seq[end], sai);
    while (sai.x[2] != 0) {
      end++;
      rldintv_t osai[6];
      rld_extend(index, &sai, osai, 0);
      sai = osai[fm6_comp(seq[end])];
    }

    // add solution
    printf("%s\t%d\t%d\n", first ? qname : "*", begin, end - begin + 1);
    first = 0;

    // prepare for next round
    if (begin == 0)
      break;
    begin = end - 1;
  }
}

int main_pingpong(int argc, char *argv[]) {
  (void)argc; // suppress unused parameter warning

  char *index_fn = argv[1];
  char *fq_fn = argv[2];

  rld_t *index = rld_restore(index_fn);

  gzFile fp = gzopen(fq_fn, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  uint8_t *s;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

    pp(index, s, l, ks->name.s);
  }
  kseq_destroy(ks);
  gzclose(fp);
  rld_destroy(index);

  return 0;
}
