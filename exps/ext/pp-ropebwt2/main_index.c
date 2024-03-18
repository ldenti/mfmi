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
  int block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES,
      so = MR_SO_RCLO;

  int64_t m = (int64_t)(.97 * 10 * 1024 * 1024 * 1024) + 1;
  int threads = 0;
  int c;
  while ((c = getopt(argc, argv, "@")) >= 0) {
    switch (c) {
    case '@':
      threads = 1;
      continue;
    default:
      return 1;
    }
  }
  if (argc - optind < 1) {
    return 1;
  }

  double t_start = realtime();

  mrope_t *mr = mr_init(max_nodes, block_len, so);
  if (threads > 0)
    mr_thr_min(mr, 10);

  int i, l;
  char *fa_path;
  gzFile fp;
  kseq_t *ks;
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  double ct, rt;
  while (optind < argc) {
    fa_path = argv[optind++];

    fp = gzopen(fa_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
      s = (uint8_t *)ks->seq.s;

      // change encoding
      for (i = 0; i < l; ++i)
        s[i] = fm6(s[i]);
      // Reverse the sequence
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        s[l - 1 - i] = s[i];
        s[i] = tmp;
      }

      // Add forward to buffer
      kputsn((char *)s, l + 1, &buf);

      // Add reverse to buffer
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        tmp = fm6_comp(tmp);
        s[l - 1 - i] = fm6_comp(s[i]);
        s[i] = tmp;
      }
      if (l & 1)
        s[i] = fm6_comp(s[i]);
      kputsn((char *)s, l + 1, &buf);

      if (buf.l >= m) {
        ct = cputime(), rt = realtime();
        mr_insert_multi(mr, buf.l, (const uint8_t *)buf.s, threads);
        fprintf(stderr,
                "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
                __func__, (long)buf.l, realtime() - rt, cputime() - ct);
        buf.l = 0;
      }
    }
    if (buf.l) {
      ct = cputime(), rt = realtime();
      mr_insert_multi(mr, buf.l, (const uint8_t *)buf.s, threads);
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
  fprintf(stderr,
          "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
          __func__, realtime() - t_start, cputime());

  mr_destroy(mr);

  return 0;
}

int main_index_v2(int argc, char *argv[]) {
  int block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES,
      so = MR_SO_RCLO;

  int64_t m = INT32_MAX;
  int threads = 0;
  int c;
  while ((c = getopt(argc, argv, "@")) >= 0) {
    switch (c) {
    case '@':
      threads = 1;
      continue;
    default:
      return 1;
    }
  }
  if (argc - optind < 1) {
    return 1;
  }

  double t_start = realtime();

  mrope_t *mr = mr_init(max_nodes, block_len, so);
  if (threads > 0)
    mr_thr_min(mr, 10);

  int i, l;
  char *fa_path;
  gzFile fp;
  kseq_t *ks;
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  kstring_t buf_rc = {0, 0, 0};
  double ct, rt;
  while (optind < argc) {
    fa_path = argv[optind++];

    fp = gzopen(fa_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
      if (buf.l + l > m) {
        ct = cputime(), rt = realtime();
        mr_insert_multi(mr, buf.l, (const uint8_t *)buf.s, threads);
        fprintf(stderr,
                "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
                __func__, (long)buf.l, realtime() - rt, cputime() - ct);
        buf.l = 0;

        ct = cputime(), rt = realtime();
        mr_insert_multi(mr, buf_rc.l, (const uint8_t *)buf_rc.s, threads);
        fprintf(stderr,
                "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
                __func__, (long)buf_rc.l, realtime() - rt, cputime() - ct);
        buf_rc.l = 0;
      }

      s = (uint8_t *)ks->seq.s;

      // change encoding
      for (i = 0; i < l; ++i)
        s[i] = fm6(s[i]);
      // Reverse the sequence
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        s[l - 1 - i] = s[i];
        s[i] = tmp;
      }

      // Add forward to buffer
      kputsn((char *)s, l + 1, &buf);

      // Add reverse to buffer
      for (i = 0; i < (l >> 1); ++i) {
        int tmp = s[l - 1 - i];
        tmp = fm6_comp(tmp);
        s[l - 1 - i] = fm6_comp(s[i]);
        s[i] = tmp;
      }
      if (l & 1)
        s[i] = fm6_comp(s[i]);
      kputsn((char *)s, l + 1, &buf_rc);
    }
    if (buf.l) {
      ct = cputime(), rt = realtime();
      mr_insert_multi(mr, buf.l, (const uint8_t *)buf.s, threads);
      fprintf(stderr,
              "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
              __func__, (long)buf.l, realtime() - rt, cputime() - ct);
      buf.l = 0;

      ct = cputime(), rt = realtime();
      mr_insert_multi(mr, buf_rc.l, (const uint8_t *)buf_rc.s, threads);
      fprintf(stderr,
              "[M::%s] inserted %ld symbols in %.3f sec; %.3f CPU sec\n",
              __func__, (long)buf_rc.l, realtime() - rt, cputime() - ct);
      buf_rc.l = 0;
    }
    kseq_destroy(ks);
    gzclose(fp);

    fprintf(stderr,
            "[M::%s] indexed %s - Total time: %.3f sec; CPU: %.3f sec\n",
            __func__, fa_path, realtime() - t_start, cputime());
  }

  free(buf.s);
  free(buf_rc.s);

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
  fprintf(stderr,
          "[M::%s] dumped index - Total time: %.3f sec; CPU: %.3f sec\n",
          __func__, realtime() - t_start, cputime());

  mr_destroy(mr);

  return 0;
}
