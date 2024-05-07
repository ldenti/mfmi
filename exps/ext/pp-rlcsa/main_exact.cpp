#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <zlib.h>

#include "fmd.h"
#include "kseq.h"
#include "rlcsa.h"

KSEQ_INIT(gzFile, gzread)

int main_exact(int argc, char **argv) {
  (void)argc; // suppress unused parameter warning

  char *index_fn = argv[1]; // index
  char *fq_fn = argv[2];    // **perfect** reads

  CSA::FMD fmd(index_fn, false);

  int errors = 0;

  gzFile fp = gzopen(fq_fn, "rb");
  kseq_t *ks = kseq_init(fp);
  int l;
  while ((l = kseq_read(ks)) >= 0) {
    auto r = fmd.fmdCount(ks->seq.s, true);
    if (r.end_offset >= 0)
      continue;
    r = fmd.fmdCount(ks->seq.s, false);
    errors += r.end_offset < 0;
  }
  printf("Errors: %d\n", errors);
  kseq_destroy(ks);
  gzclose(fp);

  return 0;
}
