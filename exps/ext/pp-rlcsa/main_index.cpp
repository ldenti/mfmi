#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "rlcsa_builder.h"

KSEQ_INIT(gzFile, gzread)

static const unsigned char rc[128] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};

uint64_t get_size_fa(const char *fa_path) {
  gzFile fp = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fp);
  uint64_t tot_l = 0;
  int l;
  while ((l = kseq_read(seq)) >= 0)
    tot_l += l + 1;
  kseq_destroy(seq);
  gzclose(fp);
  return tot_l;
}

void concatenate_fa(const char *fa_path, uint64_t size, unsigned char *data,
                    unsigned char *data_r) {
  gzFile fp = gzopen(fa_path, "r");
  kseq_t *ks = kseq_init(fp);
  int l, i;
  uint64_t curr_l = 0;
  while ((l = kseq_read(ks)) >= 0) {
    memmove(data + curr_l, ks->seq.s, l);

    // reverse
    for (i = 0; i < (l >> 1); ++i) {
      int tmp = ks->seq.s[l - 1 - i];
      tmp = rc[tmp];
      ks->seq.s[l - 1 - i] = rc[ks->seq.s[i]];
      ks->seq.s[i] = tmp;
    }
    if (l & 1)
      ks->seq.s[i] = rc[ks->seq.s[i]];
    memmove(data_r + curr_l, ks->seq.s, l);

    curr_l += l + 1;
  }
  kseq_destroy(ks);
  gzclose(fp);
}

int main_index(int argc, char **argv) {
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
  std::vector<std::string> inputs;
  while (optind < argc)
    inputs.push_back(argv[optind++]);

  std::cerr << "Initializing.." << std::endl;
  CSA::RLCSABuilder builder(block_size, sample_rate, 0, threads, NULL);

  uint64_t max_size = ((uint64_t)1 << 32) - 1;
  unsigned char *data =
      (unsigned char *)calloc(max_size, sizeof(unsigned char));
  unsigned char *data_r =
      (unsigned char *)calloc(max_size, sizeof(unsigned char));
  for (const std::string &fa_path : inputs) {
    std::cerr << "Reading " << fa_path << std::endl;
    uint64_t size = get_size_fa(fa_path.c_str());
    concatenate_fa(fa_path.c_str(), size, data, data_r);

    std::cerr << "Indexing " << fa_path << std::endl;
    CSA::RLCSA *index =
        new CSA::RLCSA(data, size, block_size, sample_rate, threads, false);
    std::cerr << "Storing " << fa_path << std::endl;
    if (index != 0 && index->isOk())
      index->writeTo(fa_path);
    delete index;
    std::cerr << "Merging " << fa_path << std::endl;
    builder.insertFromFile(fa_path, data);

    std::cerr << "Indexing " << fa_path << " (R)" << std::endl;
    index =
        new CSA::RLCSA(data_r, size, block_size, sample_rate, threads, false);
    std::cerr << "Storing " << fa_path << " (R)" << std::endl;
    char *rev = (char *)malloc(fa_path.size() + 2);
    strcpy(rev, fa_path.c_str());
    strcat(rev, ".R");
    if (index != 0 && index->isOk())
      index->writeTo(rev);
    delete index;
    std::cerr << "Merging " << fa_path << " (R)" << std::endl;
    builder.insertFromFile(rev, data_r);
  }
  free(data);

  std::cerr << "Storing full index to " << index_prefix << std::endl;
  CSA::RLCSA *rlcsa = builder.getRLCSA();
  if (rlcsa != 0 && rlcsa->isOk()) {
    std::cout << std::endl;
    rlcsa->printInfo();
    rlcsa->reportSize(true);
    rlcsa->writeTo(index_prefix);
  }
  delete rlcsa;

  return 0;
}
