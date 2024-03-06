#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

int main_index(int argc, char *argv[]) {
  (void)argc; // suppress unused parameter warning
  double t_start = realtime();

  int nt = 1;
  int reverse = 0;
  uint32_t m = (uint32_t)(0.85 * UINT32_MAX) + 1;

  int c;
  while ((c = getopt(argc, argv, "@:rh")) >= 0) {
    switch (c) {
    case '@':
      nt = atoi(optarg);
      continue;
    case 'r':
      reverse = 1;
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

  int i, l;
  char *fa_path;
  rlcsa_t *rlc = rlc_init();
  gzFile fp;
  kseq_t *ks;
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  while (optind < argc) {
    fa_path = argv[optind++];

    fp = gzopen(fa_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
      s = (uint8_t *)ks->seq.s;

      // change encoding
      for (i = 0; i < l; ++i)
        s[i] = fm6(s[i]);

      // Add forward to buffer
      kputsn((char *)ks->seq.s, ks->seq.l + 1, &buf);
      if (reverse) {
        // Add reverse to buffer
        for (i = 0; i < (l >> 1); ++i) {
          int tmp = s[l - 1 - i];
          tmp = fm6_comp(tmp);
          s[l - 1 - i] = fm6_comp(s[i]);
          s[i] = tmp;
        }
        if (l & 1)
          s[i] = fm6_comp(s[i]);
        kputsn((char *)s, ks->seq.l + 1, &buf);
      }
      if (buf.l >= m) {
        double ct = cputime(), rt = realtime();
        rlc_insert(rlc, (const uint8_t *)buf.s, (int64_t)buf.l, nt);
        fprintf(stderr,
                "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
                __func__, (long)buf.l, realtime() - rt, cputime() - ct);
        buf.l = 0;
      }
    }
    if (buf.l) {
      double ct = cputime(), rt = realtime();
      rlc_insert(rlc, (const uint8_t *)buf.s, (int64_t)buf.l, nt);
      fprintf(stderr,
              "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
              __func__, (long)buf.l, realtime() - rt, cputime() - ct);
      buf.l = 0;
    }
    kseq_destroy(ks);
    gzclose(fp);

    fprintf(stderr,
            "[M::%s] indexed %s - Total time: %.3f sec; CPU: %.3f sec\n",
            __func__, fa_path, realtime() - t_start, cputime());
  }
  free(buf.s);

  // Print BWT - debug
  // rpitr_t *it = calloc(1, sizeof(rpitr_t));
  // rope_itr_first(rlc->rope, it);
  // uint8_t *b;
  // while ((b = (uint8_t *)rope_itr_next_block(it)) != 0)
  //   rle_print(b, 1);

  rlc_dump(rlc, "-"); // TODO: add path to CLI
  rlc_destroy(rlc);

  fprintf(stderr,
          "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
          __func__, realtime() - t_start, cputime());

  return 0;
}

int main_index_single(int argc, char *argv[]) {
  (void)argc; // suppress unused parameter warning
  double t_start = realtime();

  int nt = 1;
  int reverse = 0;
  int c;
  while ((c = getopt(argc, argv, "@:rh")) >= 0) {
    switch (c) {
    case '@':
      nt = atoi(optarg);
      continue;
    case 'r':
      reverse = 1;
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

  int i, l;
  char *fa_path;
  rlcsa_t *rlc = rlc_init();
  gzFile fp;
  kseq_t *ks;
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  kstring_t buf_rc = {0, 0, 0};
  while (optind < argc) {
    fa_path = argv[optind++];

    fp = gzopen(fa_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
      s = (uint8_t *)ks->seq.s;

      // change encoding
      for (i = 0; i < l; ++i)
        s[i] = fm6(s[i]);

      // Add forward to buffer
      kputsn((char *)ks->seq.s, ks->seq.l + 1, &buf);
      if (reverse) {
        // Add reverse to buffer
        for (i = 0; i < (l >> 1); ++i) {
          int tmp = s[l - 1 - i];
          tmp = fm6_comp(tmp);
          s[l - 1 - i] = fm6_comp(s[i]);
          s[i] = tmp;
        }
        if (l & 1)
          s[i] = fm6_comp(s[i]);
        kputsn((char *)s, ks->seq.l + 1, &buf_rc);
      }
    }
    double ct = cputime(), rt = realtime();
    rlc_insert(rlc, (const uint8_t *)buf.s, (int64_t)buf.l, nt);
    fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
            __func__, (long)buf.l, realtime() - rt, cputime() - ct);
    buf.l = 0;

    ct = cputime(), rt = realtime();
    rlc_insert(rlc, (const uint8_t *)buf_rc.s, (int64_t)buf_rc.l, nt);
    fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
            __func__, (long)buf_rc.l, realtime() - rt, cputime() - ct);
    buf_rc.l = 0;

    kseq_destroy(ks);
    gzclose(fp);

    fprintf(stderr,
            "[M::%s] indexed %s - Total time: %.3f sec; CPU: %.3f sec\n",
            __func__, fa_path, realtime() - t_start, cputime());
  }
  free(buf.s);
  free(buf_rc.s);

  rlc_dump(rlc, "-");
  rlc_destroy(rlc);

  fprintf(stderr,
          "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
          __func__, realtime() - t_start, cputime());

  return 0;
}
