#include <string.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.hpp"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

int main_index(int argc, char *argv[]);
// int main_search(int argc, char *argv[]);
// int main_pingpong(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "index") == 0)
    return main_index(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "search") == 0)
  //   return main_search(argc - 1, argv + 1);
  // else if (strcmp(argv[1], "pingpong") == 0)
  //   return main_pingpong(argc - 1, argv + 1);
  else {
    gzFile fp = gzopen(argv[1], "rb");
    kseq_t *ks = kseq_init(fp);
    uint8_t *s;
    kstring_t buf = {0, 0, 0};
    int i, l;

    rlcsa_t *rlc = rlc_init();

    l = kseq_read(ks);
    printf("Seq: %s\n", ks->seq.s);
    s = (uint8_t *)ks->seq.s;
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
    kputsn((char *)s, l + 1, &buf);
    rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);

    buf.l = 0;
    l = kseq_read(ks);
    printf("Seq: %s\n", ks->seq.s);
    s = (uint8_t *)ks->seq.s;
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
    kputsn((char *)s, l + 1, &buf);
    rlc_insert(rlc, (const uint8_t *)buf.s, (uint32_t)buf.l, 1);

    rlc_destroy(rlc);

    kseq_destroy(ks);
    gzclose(fp);

    // return 1;
  }
}
