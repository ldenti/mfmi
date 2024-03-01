#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "mrope.h"
#include "rld0.h"
#include "rle.h"

#include "utils.h"

KSEQ_INIT(gzFile, gzread)

int main_index(int argc, char *argv[]) {
  // hardcoded parameters
  int block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES,
      so = MR_SO_RCLO;

  // the index
  mrope_t *mr = mr_init(max_nodes, block_len, so);

  // Parsing the input sample
  int optind = 1;
  int i, l;
  char *fa_path;
  gzFile fp;
  kseq_t *ks;
  uint8_t *s;
  while (optind < argc) {
    fa_path = argv[optind++];

    fp = gzopen(fa_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
      s = (uint8_t *)ks->seq.s;

      // change encoding
      for (i = 0; i < l; ++i)
        s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

      // Reverse the sequence
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        s[l - 1 - i] = s[i];
        s[i] = tmp;
      }

      mr_insert1(mr, s);

      // Add reverse to buffer
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
        s[l - 1 - i] = fm6_comp(s[i]);
        s[i] = tmp;
      }
      if (l & 1)
        s[i] = fm6_comp(s[i]);
      mr_insert1(mr, s);
    }
    kseq_destroy(ks);
    gzclose(fp);
  }

  // dump index to stdout
  mritr_t itr;
  const uint8_t *block;
  rld_t *e = 0;
  rlditr_t di;
  e = rld_init(6, 3);
  rld_itr_init(e, &di, 0);
  mr_itr_first(mr, &itr, 1);
  while ((block = mr_itr_next_block(&itr)) != 0) {
    const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
    while (q < end) {
      int c = 0;
      int64_t l;
      rle_dec1(q, c, l);
      rld_enc(e, &di, l, c);
    }
  }
  rld_enc_finish(e, &di);
  rld_dump(e, "-");

  mr_destroy(mr);

  return 0;
}
