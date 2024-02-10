#include "rlcsa.h"

rlcsa_t *rlc_init() {
  rlcsa_t *rlc;
  rlc = calloc(0, sizeof(rlcsa_t));
  rlc->cnts = calloc(5, sizeof(uint64_t));
  rlc->bits = calloc(5, sizeof(uint8_t));
  return rlc;
}

int rlc_insert(rlcsa_t *rlc, uint8_t *seq, int n) {
  int i;
  for (i = 0; i < n; ++i)
    ++rlc->cnts[seq[i] - 1];

  sa_t *sa = 0;
  // sa = simpleSuffixSort(data, bytes, this->number_of_sequences, threads);
  return 0;
}

void rlc_destroy(rlcsa_t *rlc) {
  free(rlc->cnts);
  free(rlc->bits);
  free(rlc);
}

// short_pair *simpleSuffixSort(const usint *sequence, uint n, uint threads) {
//   if (sequence == 0 || n == 0) {
//     return 0;
//   }

//   skew_pair *pairs =
//       (skew_pair
//            *)new short_pair[n * sizeof(skew_pair) / sizeof(short_pair) + 1];
//   uint *keys = new uint[n]; // In text order.
//   std::vector<ss_range> unsorted;
//   threads = std::max(threads, (uint)1);
//   omp_set_num_threads(threads);

// // Initialize pairs.
// #pragma omp parallel for schedule(static)
//   for (uint i = 0; i < n; i++) {
//     pairs[i].first = i;
//     pairs[i].second = sequence[i];
//   }

//   // Sort according to first character.
//   parallelSort(pairs, pairs + n, skew_comparator);
//   unsorted.push_back(ss_range(0, n - 1));
//   uint total = setRanks(pairs, keys, n, unsorted, threads, 1);

//   if (sizeof(usint) < 2 * sizeof(uint)) {
//     return prefixDoubling(packPairs(pairs, n), keys, unsorted, n, threads,
//                           total, 1);
//   } else {
//     return prefixTripling(pairs, keys, unsorted, n, threads, total, 1);
//   }
// }

// short_pair *simpleSuffixSort(const uchar *sequence, uint n, uint sequences,
//                              uint threads) {
//   if (sequence == 0 || n == 0) {
//     return 0;
//   }

//   short_pair *pairs = new short_pair[n]; // In sorted order.
//   uint *keys = new uint[n];              // In text order.
//   std::vector<ss_range> unsorted;
//   threads = std::max(threads, (uint)1);
//   omp_set_num_threads(threads);

//   // Remap alphabet.
//   uint alphabet[CHARS];
//   for (uint c = 0; c < CHARS; c++) {
//     alphabet[c] = 0;
//   }
//   for (uint i = 0; i < n; i++) {
//     alphabet[sequence[i]]++;
//   }
//   uint alphabet_size = sequences;
//   for (uint c = 1; c < CHARS; c++) {
//     uint temp = alphabet_size;
//     if (alphabet[c] > 0) {
//       alphabet_size++;
//     }
//     alphabet[c] = temp;
//   }

//   // Determine pack factor.
//   uint limit = std::numeric_limits<uint>::max() / alphabet_size;
//   uint h = 1, pack_multiplier = 1;
//   if (alphabet_size > 1) {
//     while (pack_multiplier * alphabet_size <= limit) {
//       h++;
//       pack_multiplier *= alphabet_size;
//     }
//   }

//   // Initialize pairs.
//   uint zeros = 0, value = 0;
//   for (uint i = 0; i < h; i++) {
//     value *= alphabet_size;
//     if (sequence[i] == 0) {
//       value += zeros;
//       zeros++;
//     } else {
//       value += alphabet[sequence[i]];
//     }
//   }
//   for (uint i = 0; i < n - h; i++) {
//     pairs[i].first = i;
//     pairs[i].second = value;
//     value = (value % pack_multiplier) * alphabet_size;
//     if (sequence[i + h] == 0) {
//       value += zeros;
//       zeros++;
//     } else {
//       value += alphabet[sequence[i + h]];
//     }
//   }
//   for (uint i = n - h; i < n; i++) {
//     pairs[i].first = i;
//     pairs[i].second = value;
//     value = (value % pack_multiplier) * alphabet_size;
//   }

//   uint total = initialSort(pairs, keys, unsorted, n, threads, h);
//   return prefixDoubling(pairs, keys, unsorted, n, threads, total, h);
// }