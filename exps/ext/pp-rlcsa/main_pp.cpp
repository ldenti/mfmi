#include <iostream>
#include <zlib.h>

#include "kseq.h"

#include "fmd.h"

KSEQ_INIT(gzFile, gzread)

void pp(const CSA::FMD *index, char *seq, const int64_t l, const char *qname) {
  bool first = true;

  int begin = l - 1;
  uint8_t c;
  CSA::FMDPosition sai;
  while (begin >= 0) {
    c = seq[begin];
    sai = index->getCharPosition(c);

    // Backward search. Find a mismatching sequence. Stop at first mismatch.
    while (!sai.isEmpty() && begin > 0) {
      --begin;
      c = seq[begin];
      sai = index->extend(sai, c, true);
    }
    // prefix is a match
    if (begin == 0 && !sai.isEmpty())
      break;

    // Forward search
    int end = begin;
    // std::cerr << "F from " << end << std::endl;
    c = seq[end];
    sai = index->getCharPosition(c);
    while (!sai.isEmpty()) {
      ++end;
      c = seq[end];
      sai = index->extend(sai, c, false);
    }

    // print solution
    std::cout << (first ? qname : "*") << "\t" << begin << "\t"
              << end << "\t" << end - begin + 1 << std::endl;
    first = false;

    // prepare for next round
    if (begin == 0)
      break;
    begin = end - 1;
  }
}

int main_pingpong(int argc, char *argv[]) {
  (void)argc; // suppress unused parameter warning

  char *index_prefix = argv[1];
  char *fq_fn = argv[2];

  const CSA::FMD *index = new CSA::FMD(index_prefix, false);

  gzFile fp = gzopen(fq_fn, "rb");
  kseq_t *ks = kseq_init(fp);
  int l;
  while ((l = kseq_read(ks)) >= 0) {
    pp(index, ks->seq.s, l, ks->name.s);
  }
  kseq_destroy(ks);
  gzclose(fp);

  delete index;

  return 0;
}
