#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "rlcsa_builder.h"

#include "utils.h"

KSEQ_INIT(gzFile, gzread)

int main_index(int argc, char *argv[]) {
  uint sample_rate = 0; // (1 << 31);
  int block_size = 32;
  int threads = 1;
  std::string index_prefix = "RLCSA";

  int c;
  while ((c = getopt(argc, argv, "i:@:")) >= 0) {
    switch (c) {
    case 'i':
      index_prefix = optarg;
      continue;
    case '@':
      threads = atoi(optarg);
      continue;
    default:
      return 1;
    }
  }
  if (argc - optind < 1) {
    return 1;
  }

  CSA::RLCSABuilder builder(block_size, sample_rate, 0, threads, NULL);

  gzFile fp;
  kseq_t *ks;
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  kstring_t buf_rc = {0, 0, 0};
  int i, l;
  char *fa_path;
  CSA::RLCSA *index;
  double rt;
  while (optind < argc) {
    fa_path = argv[optind++];
    fp = gzopen(fa_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
      if (buf.l + l > INT32_MAX) {
        rt = realtime();
        index = new CSA::RLCSA((CSA::uchar *)buf.s, buf.l, block_size,
                               sample_rate, threads, false);
        index->writeTo(fa_path); // TODO: store to tmp dir
        delete index;
        builder.insertFromFile(fa_path, (CSA::uchar *)buf.s);
        fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec\n", __func__,
                (long)buf.l, realtime() - rt);
        buf.l = 0;

        rt = realtime();
        index = new CSA::RLCSA((CSA::uchar *)buf_rc.s, buf_rc.l, block_size,
                               sample_rate, threads, false);
        index->writeTo(fa_path); // TODO: store to tmp dir
        delete index;
        builder.insertFromFile(fa_path, (CSA::uchar *)buf_rc.s);
        fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec\n", __func__,
                (long)buf_rc.l, realtime() - rt);
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
      rt = realtime();
      index = new CSA::RLCSA((CSA::uchar *)buf.s, buf.l, block_size,
                             sample_rate, threads, false);
      index->writeTo(fa_path); // TODO: store to tmp dir
      delete index;
      builder.insertFromFile(fa_path, (CSA::uchar *)buf.s);
      fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec\n", __func__,
              (long)buf.l, realtime() - rt);
      buf.l = 0;

      rt = realtime();
      index = new CSA::RLCSA((CSA::uchar *)buf_rc.s, buf_rc.l, block_size,
                             sample_rate, threads, false);
      index->writeTo(fa_path); // TODO: store to tmp dir
      delete index;
      builder.insertFromFile(fa_path, (CSA::uchar *)buf_rc.s);
      fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec\n", __func__,
              (long)buf_rc.l, realtime() - rt);
      buf_rc.l = 0;
    }
    kseq_destroy(ks);
    gzclose(fp);
  }
  free(buf.s);
  free(buf_rc.s);

  rt = realtime();
  index = builder.getRLCSA();
  if (index != 0 && index->isOk()) {
    // std::cout << std::endl;
    // rlcsa->printInfo();
    // rlcsa->reportSize(true);
    index->writeTo(index_prefix);
  }
  fprintf(stderr, "[M::%s] dumped index in %.3f sec\n", __func__,
          realtime() - rt);
  delete index;

  return 0;
}
