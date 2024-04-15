#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "utils.h"
#include "rld0.h"

KSEQ_INIT(gzFile, gzread)

int main_index(int argc, char *argv[]) {
  int nt = 1;
  // uint32_t m = (uint32_t)(INT32_MAX); // (uint32_t)(0.85 * UINT32_MAX) + 1;
  char otype = 'r'; // rope
  int c;
  while ((c = getopt(argc, argv, "@:dh")) >= 0) {
    switch (c) {
    case '@':
      nt = atoi(optarg);
      continue;
    case 'd':
      otype = 'd'; // delta
      continue;
    case 'h':
      printf("HELP\n");
      return 0;
    default:
      printf("HELP\n");
      return 1;
    }
  }
  if (argc - optind < 1) {
    printf("HELP\n");
    return 1;
  }

  rlcsa_t *rlc = rlc_init();

  gzFile fp;
  kseq_t *ks;
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  kstring_t buf_rc = {0, 0, 0};
  int i, l;
  char *fa_path;
  double ct, rt;
  while (optind < argc) {
    fa_path = argv[optind++];
    fp = gzopen(fa_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
      if (buf.l + l > INT32_MAX) {
        ct = cputime(), rt = realtime();
        rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, nt);
        fprintf(stderr,
                "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
                __func__, (long)buf.l, realtime() - rt, cputime() - ct);
        buf.l = 0;
        ct = cputime(), rt = realtime();
        rlc_insert(rlc, (const uint8_t *)buf_rc.s, (uint32_t)buf_rc.l, nt);
        fprintf(stderr,
                "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
                __func__, (long)buf_rc.l, realtime() - rt, cputime() - ct);
        buf_rc.l = 0;
      }

      s = (uint8_t *)ks->seq.s;

      // change encoding
      for (i = 0; i < l; ++i)
        s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

      // Add forward to buffer
      kputsn((char *)s, l + 1, &buf);

      // Add reverse to buffer
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
        s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
        s[i] = tmp;
      }
      if (l & 1)
        s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
      kputsn((char *)s, ks->seq.l + 1, &buf_rc);
    }
    if (buf.l) {
      ct = cputime(), rt = realtime();
      rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, nt);
      fprintf(stderr,
              "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
              __func__, (long)buf.l, realtime() - rt, cputime() - ct);
      buf.l = 0;
      ct = cputime(), rt = realtime();
      rlc_insert(rlc, (const uint8_t *)buf_rc.s, (uint32_t)buf_rc.l, nt);
      fprintf(stderr,
              "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
              __func__, (long)buf_rc.l, realtime() - rt, cputime() - ct);
      buf_rc.l = 0;
    }
    kseq_destroy(ks);
    gzclose(fp);
  }
  free(buf.s);
  free(buf_rc.s);

  // Print BWT - debug
  // rpitr_t *it = calloc(1, sizeof(rpitr_t));
  // rope_itr_first(rlc->rope, it);
  // uint8_t *b;
  // while ((b = (uint8_t *)rope_itr_next_block(it)) != 0)
  //   rle_print(b, 1);

  ct = cputime(), rt = realtime();
  if (otype == 'r')
    rlc_dump(rlc, "-"); // TODO: add path to CLI
  else {
    rpitr_t itr;
    const uint8_t *block;
    rld_t *e = 0;
    rlditr_t di;
    e = rld_init(6, 3);
    rld_itr_init(e, &di, 0);
    rope_itr_first(rlc->rope, &itr);
    while ((block = rope_itr_next_block(&itr)) != 0) {
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
  }
  rlc_destroy(rlc);
  fprintf(stderr, "[M::%s] dumped index in  %.3f sec; %.3f CPU sec\n", __func__,
          realtime() - rt, cputime() - ct);

  return 0;
}
