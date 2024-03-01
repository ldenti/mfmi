#include <iostream>
#include <zlib.h>

#include "kseq.h"

#include "fmd_simple.hpp"

KSEQ_INIT(gzFile, gzread)

std::string interval2str(const FMDPosition &i) {
  return std::to_string(i.forward_start) + "," +
         std::to_string(i.reverse_start) + "," +
         std::to_string(i.end_offset + 1);
}

void pp(const FMD *index, uint8_t *seq, const int64_t l, const char *qname) {
  bool first = true;

  int begin = l - 1;
  uint8_t c;
  FMDPosition sai;
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
    c = seq[end];
    sai = index->getCharPosition(c);
    while (!sai.isEmpty()) {
      ++end;
      c = seq[end];
      sai = index->extend(sai, c, false);
    }

    // add solution
    std::cout << (first ? qname : "*") << "\t" << begin << "\t"
              << end - begin + 1 << std::endl;
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

  const FMD *index = new FMD(index_prefix, false);

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
  delete index;

  return 0;
}
