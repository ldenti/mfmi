#ifndef RLCSA_H_
#define RLCSA_H_

#include "rle.h"
#include <stdlib.h>

typedef struct {
  uint a;
  uint b;
} sa_t;

typedef struct {
  uint64_t *cnts;
  uint8_t *bits;
} rlcsa_t;

rlcsa_t *rlc_init();

int rlc_insert(rlcsa_t *rlc, uint8_t *seq, int n);

void rlc_destroy(rlcsa_t *rlc);

#endif

// void RLCSA::buildRLCSA(uchar *data, usint *ranks, usint bytes, usint
// block_size,
//                        usint threads) {

//   threads = std::max(threads, (usint)1);
//   omp_set_num_threads(threads);

//   // Determine the number of sequences and mark their end points.
//   DeltaEncoder endings(RLCSA::ENDPOINT_BLOCK_SIZE);

//   this->number_of_sequences = 0;
//   usint marker = 0;
//   usint padding = 0, chars_encountered = 0;

//   for (usint i = 0; i < bytes; i++) {
//     if (data[i] == 0) {
//       if (i == marker) {
//         break;
//       } // Empty sequence.
//       this->number_of_sequences++;
//       marker = i + 1;
//       usint pos = chars_encountered + padding - 1;
//       endings.setBit(pos);
//       padding = ((pos + 1) / 1) * 1 - chars_encountered;
//     } else {
//       ++chars_encountered;
//     }
//   }

//   if (this->number_of_sequences == 0 || marker != bytes) {
//     std::cerr << "RLCSA: Collection must consist of 0-terminated nonempty "
//                  "sequences !"
//               << std::endl;
//     return;
//   }
//   this->end_points = new DeltaVector(endings, chars_encountered + padding);

//   // Build character tables etc.
//   usint distribution[CHARS];
//   for (usint c = 0; c < CHARS; c++) {
//     distribution[c] = 0;
//   }
//   for (usint i = 0; i < bytes; i++) {
//     distribution[(usint)data[i]]++;
//   }
//   distribution[0] = 0; // \0 is an end marker
//   this->alphabet = new Alphabet(distribution);
//   this->data_size = this->alphabet->getDataSize();

//   // Build suffix array.
//   short_pair *sa = 0;
//   if (ranks == 0) {
//     sa = simpleSuffixSort(data, bytes, this->number_of_sequences, threads);
//   } else {
//     sa = simpleSuffixSort(ranks, bytes, threads);
//   }

// // Build Psi.
// #pragma omp parallel for schedule(static)
//   for (usint i = 0; i < bytes; i++) {
//     sa[i].first = sa[(sa[i].first + 1) % bytes].second;
//   }

// // Build RLCSA.
// #pragma omp parallel for schedule(dynamic, 1)
//   for (usint c = 0; c < CHARS; c++) {
//     if (!(this->alphabet->hasChar(c))) {
//       this->array[c] = 0;
//       continue;
//     }

//     short_pair *curr =
//         sa + this->alphabet->cumulative(c) + this->number_of_sequences;
//     short_pair *limit = curr + this->alphabet->countOf(c);

//     sdsl::rle_vector_builder<64> encoder(this->data_size +
//                                          this->number_of_sequences);
//     pair_type run((*curr).first, 1);
//     ++curr;

//     for (; curr < limit; ++curr) {
//       if ((*curr).first == run.first + run.second) {
//         run.second++;
//       } else {
//         encoder.set(run.first, run.second);
//         run = pair_type((*curr).first, 1);
//       }
//     }
//     encoder.set(run.first, run.second);
//     // encoder.flush();

//     this->array[c] = new sdsl::rle_vector<64>(encoder);
//   }
//   delete[] sa;

//   this->ok = true;
// }