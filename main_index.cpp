#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.hpp"
#include "rld0.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

int main_index(int argc, char *argv[]) {
  int nt = 1;
  int c;
  while ((c = getopt(argc, argv, "@:dh")) >= 0) {
    switch (c) {
    case '@':
      nt = atoi(optarg);
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

  // Build BWT and store as rld0
  ct = cputime(), rt = realtime();
  // rlc_print_bwt(rlc);
  rlc_dump(rlc, "-");

  rlc_destroy(rlc);

  fprintf(stderr, "[M::%s] dumped index in  %.3f sec; %.3f CPU sec\n", __func__,
          realtime() - rt, cputime() - ct);

  return 0;
}
