#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "rld0.h"

#include "utils.h"

KSEQ_INIT(gzFile, gzread)

int main_exact(int argc, char *argv[]) {
  (void)argc; // suppress unused parameter warning

  char *index_fn = argv[1]; // index
  char *fq_fn = argv[2];    // **perfect** reads

  rld_t *index = rld_restore(index_fn);

  rldintv_t sai;
  for (int c = 0; c < 6; ++c) {
    fm6_set_intv(index, c, sai);
    printf("%d: [%ld, %ld, %ld]\n", c, sai.x[0], sai.x[1], sai.x[2]);
  }

  int errors = 0;
  gzFile fp = gzopen(fq_fn, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  uint8_t *s;
  rldintv_t osai[6];
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

    // init search
    i = l - 1;
    fm6_set_intv(index, s[i], sai);
    --i;

    // backward extensions
    for (; i >= 0; --i) {
      rld_extend(index, &sai, osai, 1);
      sai = osai[s[i]];
      if (sai.x[2] <= 0) {
        errors += 1;
        break;
      }
    }
  }
  printf("Errors: %d\n", errors);
  kseq_destroy(ks);
  gzclose(fp);
  rld_destroy(index);

  return 0;
}
