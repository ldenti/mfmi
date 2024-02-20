#include "rlcsa.h"

int sat_comparator(const sa_t *a, const sa_t *b) { return a->b - b->b; }

int skew_comparator(const skew_pair *a, const skew_pair *b) {
  return a->b - b->b;
}

uint32_t setRanks(skew_pair *pairs, uint32_t *keys, uint64_t n,
                  ss_ranges *unsorted) { // uint8_t threads, uint32_t chunk) {
  ss_ranges buffer;
  kv_init(buffer);
  uint32_t total = 0;
  ss_range a;

  // TODO: parallelize
  // std::vector<ss_range> buffers[threads];
  // uint subtotals[threads];
  // for (uint i = 0; i < threads; i++) {
  //   subtotals[i] = 0;
  // }

  // #pragma omp parallel for schedule(dynamic, chunk)
  for (uint32_t i = 0; i < kv_size(*unsorted); ++i) {
    // uint thread = omp_get_thread_num();
    skew_pair *prev = pairs + kv_A(*unsorted, i).a;
    keys[prev->a] = kv_A(*unsorted, i).a;
    skew_pair *limit = pairs + kv_A(*unsorted, i).b;
    while (prev < limit) {
      skew_pair *curr = prev + 1;
      if (curr->b != prev->b) {
        ++prev;
        keys[prev->a] = prev - pairs;
        continue;
      }
      keys[curr->a] = prev - pairs;
      for (++curr; curr <= limit && curr->b == prev->b; ++curr)
        keys[curr->a] = prev - pairs;
      a = (ss_range){prev - pairs, (curr - 1) - pairs};
      kv_push(ss_range, buffer, a);
      total += curr - prev;
      // buffers[thread].push_back(ss_range(prev - pairs, (curr - 1) - pairs));
      // subtotals[thread] += curr - prev;
      prev = curr;
      if (prev <= limit)
        keys[prev->a] = prev - pairs;
    }
  }
  unsorted->n = 0;
  for (uint i = 0; i < kv_size(buffer); ++i)
    kv_push(ss_range, *unsorted, kv_A(buffer, i));
  // for (uint i = 0; i < threads; i++) {
  //   unsorted.insert(unsorted.end(), buffers[i].begin(), buffers[i].end());
  //   total += subtotals[i];
  // }
  kv_destroy(buffer);
  return total;
}

sa_t *prefixDoubling(sa_t *pairs, uint32_t *keys, ss_ranges *unsorted, uint n,
                     uint threads, uint total, uint h) {
  // Double prefix length until sorted
  // TODO: parallelize
  while (total > 0) {
    // uint chunk = std::max((size_t)1, unsorted.size() / (threads * threads));
    // #pragma omp parallel for schedule(dynamic, chunk)
    for (uint32_t i = 0; i < kv_size(*unsorted); ++i) {
      // Set sort keys for the current range.
      for (sa_t *curr = pairs + kv_A(*unsorted, i).a;
           curr <= pairs + kv_A(*unsorted, i).b; ++curr)
        curr->b = keys[curr->a + h];
      // TODO: qsort vs ksort
      qsort(pairs + kv_A(*unsorted, i).a,
            kv_A(*unsorted, i).b - kv_A(*unsorted, i).a + 1, sizeof(sa_t),
            sat_comparator); // TODO: change name of comparator
    }
    total = setRanks(pairs, keys, n, unsorted); //, threads, chunk);
    h *= 2;
    fprintf(stderr, "Sorted with %d, unsorted total = %d (%ld ranges)\n", h,
            total, kv_size(*unsorted));
  }
  // TODO: parallelize
  // #pragma omp parallel for schedule(static)
  for (uint i = 0; i < n; i++)
    pairs[i].b = keys[i];
  return pairs;
}

sa_t *simpleSuffixSort(const uint8_t *sequence, uint32_t n, uint8_t threads) {
  if (sequence == 0 || n == 0)
    return 0;
  ss_ranges unsorted;
  kv_init(unsorted);
  skew_pair *pairs = (skew_pair *)calloc(n + 1, sizeof(skew_pair));
  uint32_t *keys =
      (uint32_t *)calloc(n + 1, sizeof(uint32_t)); // In text order.
  // Initialize pairs.
  // TODO: parallelize
  // #pragma omp parallel for schedule(static)
  for (uint32_t i = 0; i < n; ++i) {
    pairs[i].a = i;
    pairs[i].b = sequence[i];
  }
  // Sort according to first character
  // TODO: ksort vs qsort - ks_mergesort(skewpair, kv_size(pairs), pairs, 0);
  qsort(pairs, n, sizeof(skew_pair), skew_comparator);
  ss_range a = (ss_range){0, n - 1};
  kv_push(ss_range, unsorted, a);
  uint32_t total = setRanks(pairs, keys, n, &unsorted); // threads, 1);
  sa_t *sa = prefixDoubling(pairs, keys, &unsorted, n, threads, total, 1);
  // TODO: doubling vs tripling?
  free(keys);
  kv_destroy(unsorted);
  return sa;
}

rlcsa_t *rlc_init(uint64_t size, uint32_t sd) {
  rlcsa_t *rlc;
  rlc = calloc(1, sizeof(rlcsa_t));
  rlc->size = size;
  rlc->cnts = calloc(6, sizeof(uint64_t));
  rlc->C = calloc(6, sizeof(uint64_t));
  for (int c = 0; c < 6; ++c)
    rlc->bits[c] = rle_init(size, sd);
  for (int c = 0; c < 6; ++c)
    rlc->end_cnts[c] = calloc(6, sizeof(int64_t));
  return rlc;
}

void rlc_destroy(rlcsa_t *rlc) {
  free(rlc->C);
  free(rlc->cnts);
  for (int c = 0; c < 6; ++c) {
    rle_destroy(rlc->bits[c]);
    free(rlc->end_cnts[c]);
  }
  free(rlc);
}

sa_t rlc_init_interval(rlcsa_t *rlc, uint8_t c) {
  return (sa_t){rlc->C[c], rlc->C[c] + rlc->cnts[c] - 1};
}

sa_t LF(rlcsa_t *rlc, sa_t range, uint8_t c) {
  int64_t cx[6] = {0, 0, 0, 0, 0, 0};
  rle_rank1a(rlc->bits[c], range.a, cx);
  int a = rlc->C[c] + cx[1];
  for (int i = 0; i < 6; ++i)
    cx[i] = 0;
  rle_rank1a(rlc->bits[c], range.b + 1, cx);
  int b = rlc->C[c] + cx[1] - 1;
  return (sa_t){a, b};
}

int rlc_insert(rlcsa_t *rlc, const uint8_t *sequence, uint n) {
  int i, c;

  // Build C
  for (i = 0; i < n; ++i)
    ++rlc->cnts[sequence[i]];

  rlc->C[0] = 0;
  for (c = 1; c < 6; ++c) // FIXME: hardcoded
    rlc->C[c] = rlc->C[c - 1] + rlc->cnts[c - 1];

  // Build SA
  sa_t *sa = simpleSuffixSort(sequence, n, 1);

  // Build Psi
  // TODO: #pragma omp parallel for schedule(static)
  for (i = 0; i < n; ++i)
    sa[i].a = sa[(sa[i].a + 1) % n].b;

  // Build RLCSA
  // TODO: #pragma omp parallel for schedule(dynamic, 1)
  for (c = 0; c < 6; ++c) { // TODO: hardcoded
    if (rlc->cnts[c] == 0)
      continue;

    sa_t *curr = sa + rlc->C[c];
    sa_t *limit = curr + rlc->cnts[c];
    uint64_t last_e = 0;
    uint64_t curr_s = curr->a;
    uint64_t curr_l = 1;
    ++curr;

    int64_t cnt[6];
    for (; curr < limit; ++curr) {
      if (curr->a == curr_s + curr_l)
        ++curr_l;
      else {
        if ((last_e == 0 && curr_s > 0) || last_e != 0) {
          rle_insert(rlc->bits[c], last_e, 0, curr_s - last_e, cnt,
                     rlc->end_cnts[c]);
          rlc->end_cnts[c][0] += curr_s - last_e;
        }
        rle_insert(rlc->bits[c], curr_s, 1, curr_l, cnt, rlc->end_cnts[c]);
        rlc->end_cnts[c][1] += curr_l;
        last_e = curr_s + curr_l;
        curr_s = curr->a;
        curr_l = 1;
      }
    }
    if ((last_e == 0 && curr_s > 0) || last_e != 0) {
      rle_insert(rlc->bits[c], last_e, 0, curr_s - last_e, cnt,
                 rlc->end_cnts[c]);
      rlc->end_cnts[c][0] += curr_s - last_e;
    }
    rle_insert(rlc->bits[c], curr_s, 1, curr_l, cnt, rlc->end_cnts[c]);
    rlc->end_cnts[c][1] += curr_l;
    curr_s += curr_l;
    if (n - curr_s != 0)
      rle_insert(rlc->bits[c], curr_s, 0, n - curr_s, cnt, rlc->end_cnts[c]);
    // rle_print(rlc->bits[c], 0);
  }
  free(sa);

  return 0;
}
