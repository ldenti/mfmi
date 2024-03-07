#include "rlcsa.h"

KSORT_INIT(pair, pair_t, pair_lt)

KSORT_INIT_GENERIC(int64_t)

rlcsa_t *rlc_init() {
  rlcsa_t *rlc;
  rlc = calloc(1, sizeof(rlcsa_t));
  rlc->cnts = calloc(6, sizeof(int64_t));
  rlc->C = calloc(6, sizeof(int64_t));
  rlc->rope = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN);
  return rlc;
}

void rlc_destroy(rlcsa_t *rlc) {
  free(rlc->cnts);
  free(rlc->C);
  rope_destroy(rlc->rope);
  free(rlc);
}

void rlc_insert(rlcsa_t *rlc, const uint8_t *sequence, uint32_t n, int nt) {
  if (rlc->l == 0) {
    rlc_build(rlc, sequence, n, nt);
    return;
  }
  rlcsa_t *rlc2 = rlc_init();
  rlc_build(rlc2, sequence, n, nt);
  rlc_merge(rlc, rlc2, sequence, nt);
}

int rlc_dump(const rlcsa_t *rlc, const char *fn) {
  FILE *fp = strcmp(fn, "-") ? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
  if (fp == 0)
    return -1;
  fwrite("RLC\3", 1, 4, fp); // write magic
  fwrite(&rlc->l, 8, 1, fp);
  fwrite(rlc->C, 8, 6, fp);
  fwrite(rlc->cnts, 8, 6, fp);
  rope_dump(rlc->rope, fp);

  // DEBUG
  // rpitr_t *it = calloc(1, sizeof(rpitr_t));
  // rope_itr_first(rlc->rope, it);
  // uint8_t *b1;
  // while ((b1 = (uint8_t *)rope_itr_next_block(it)) != 0)
  //   rle_print(b1, 1);

  fclose(fp);
  return 1;
}

rlcsa_t *rlc_restore(const char *fn) {
  rlcsa_t *rlc = rlc_init();
  FILE *fp;
  if ((fp = fopen(fn, "rb")) == 0)
    return rlc; // TODO: fail
  char magic[4];
  fread(magic, 1, 4, fp);
  fread(&rlc->l, 8, 1, fp);
  fread(rlc->C, 8, 6, fp);
  fread(rlc->cnts, 8, 6, fp);
  rlc->rope = rope_restore(fp);
  return rlc;
}

/**************************/
/** *** CONSTRUCTION *** **/
/**************************/
uint32_t setRanks(skew_pair *pairs, uint32_t *keys, ss_ranges *unsorted,
                  uint nt, uint chunk) {
  uint32_t i, p;
  ss_ranges *buffers = calloc(nt, sizeof(ss_ranges));
  for (i = 0; i < nt; ++i)
    kv_init(buffers[i]);
  uint32_t *subtotals = calloc(nt, sizeof(uint64_t));

#pragma omp parallel for num_threads(nt) schedule(dynamic, chunk)
  for (i = 0; i < kv_size(*unsorted); ++i) {
    int thread = omp_get_thread_num();
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
      ss_range a = (ss_range){prev - pairs, (curr - 1) - pairs};
      kv_push(ss_range, buffers[thread], a);
      subtotals[thread] += curr - prev;
      prev = curr;
      if (prev <= limit)
        keys[prev->a] = prev - pairs;
    }
  }
  unsorted->n = 0;
  uint32_t total = 0;
  for (i = 0; i < nt; ++i) {
    for (p = 0; p < kv_size(buffers[i]); ++p)
      kv_push(ss_range, *unsorted, kv_A(buffers[i], p));
    kv_destroy(buffers[i]);
    total += subtotals[i];
  }
  free(buffers);
  free(subtotals);
  return total;
}

sa_t *prefixDoubling(sa_t *pairs, uint32_t *keys, ss_ranges *unsorted,
                     uint32_t n, uint32_t total, int h, int nt) {
  // Double prefix length until sorted
  uint32_t i;
  while (total > 0) {
    uint chunk = max(1, kv_size(*unsorted) / (nt * nt));
#pragma omp parallel for num_threads(nt) schedule(dynamic, chunk)
    for (i = 0; i < kv_size(*unsorted); ++i) {
      // Set sort keys for the current range.
      for (sa_t *curr = pairs + kv_A(*unsorted, i).a;
           curr <= pairs + kv_A(*unsorted, i).b; ++curr)
        curr->b = keys[curr->a + h];
      ks_mergesort(pair, kv_A(*unsorted, i).b - kv_A(*unsorted, i).a + 1,
                   pairs + kv_A(*unsorted, i).a, 0);
    }
    total = setRanks(pairs, keys, unsorted, nt, chunk);
    h *= 2;
    // fprintf(stderr, "Sorted with %d, unsorted total = %ld (%ld ranges)\n", h,
    // total, kv_size(*unsorted));
  }

#pragma omp parallel for num_threads(nt) schedule(static)
  for (i = 0; i < n; i++)
    pairs[i].b = keys[i];
  return pairs;
}

sa_t *simpleSuffixSort(const uint8_t *sequence, uint32_t n, uint32_t nsep,
                       int nt) {
  if (sequence == 0 || n == 0)
    // TODO: fail more gracefully
    exit(1);

  double t = realtime();
  ss_ranges unsorted;
  kv_init(unsorted);
  skew_pair *pairs = (skew_pair *)calloc(n + 1, sizeof(skew_pair));
  uint32_t *keys =
      (uint32_t *)calloc(n + 1, sizeof(uint32_t)); // In text order.

  // Determine pack factor
  uint32_t alphabet_size = nsep + 5;           // FIXME: hardcoded
  uint32_t limit = UINT32_MAX / alphabet_size; // FIXME
  uint32_t h = 1, pack_multiplier = 1;
  while (pack_multiplier * alphabet_size <= limit) {
    ++h;
    pack_multiplier *= alphabet_size;
  }
  // FIXME: uncomment this for very short strings
  // h = h > n ? 1 : h;

  // Initialize pairs
  uint32_t zeros = 0, value = 0, i;
  for (i = 0; i < h; ++i) {
    value *= alphabet_size;
    if (sequence[i] == 0) {
      value += zeros;
      ++zeros;
    } else {
      value += nsep + sequence[i] - 1;
    }
  }
  for (i = 0; i < n - h; i++) {
    pairs[i].a = i;
    pairs[i].b = value;
    value = (value % pack_multiplier) * alphabet_size;
    if (sequence[i + h] == 0) {
      value += zeros;
      ++zeros;
    } else {
      value += nsep + sequence[i + h] - 1;
    }
  }
  for (i = n - h; i < n; i++) {
    pairs[i].a = i;
    pairs[i].b = value;
    value = (value % pack_multiplier) * alphabet_size;
  }

  // Sort according to first character
  t = realtime();
  ks_mergesort(pair, n, pairs, 0);
  fprintf(stderr, "[M::%s] Sorted pairs in %.3f sec\n", __func__,
          realtime() - t);
  t = realtime();
  ss_range a = (ss_range){0, n - 1};
  kv_push(ss_range, unsorted, a);
  uint32_t total = setRanks(pairs, keys, &unsorted, nt, 1);
  fprintf(stderr, "[M::%s] Set ranks in %.3f sec\n", __func__, realtime() - t);

  // TODO: doubling vs tripling?

  t = realtime();
  sa_t *sa = prefixDoubling(pairs, keys, &unsorted, n, total, h, nt);
  fprintf(stderr, "[M::%s] Prefix doubling in %.3f sec\n", __func__,
          realtime() - t);
  free(keys);
  kv_destroy(unsorted);
  return sa;
}

void rlc_build(rlcsa_t *rlc, const uint8_t *sequence, uint32_t n, int nt) {
  double t;
  uint32_t i, rl, p;
  uint8_t c, lc;
  rlc->l = (int64_t)n;

  // Build C
  for (i = 0; i < n; ++i)
    ++rlc->cnts[sequence[i]];
  rlc->C[0] = 0;
  for (c = 1; c < 6; ++c) // FIXME: hardcoded
    rlc->C[c] = rlc->C[c - 1] + rlc->cnts[c - 1];

  // Build SA
  t = realtime();
  sa_t *sa = simpleSuffixSort(sequence, n, (uint32_t)rlc->cnts[0], nt);
  fprintf(stderr, "[M::%s] Built SA in %.3f sec\n", __func__, realtime() - t);

  // Build rope by inserting runs
  t = realtime();
  lc = sequence[sa[0].a == 0 ? n : sa[0].a - 1]; // current run symbol
  rl = 1;                                        // current run length
  p = 0;                                         // insertion position
  for (i = 1; i < n; ++i) {
    c = sequence[sa[i].a == 0 ? n : sa[i].a - 1];
    if (c == lc) {
      ++rl;
    } else {
      rope_insert_run(rlc->rope, p, lc, rl, 0);
      p += rl;
      rl = 1;
    }
    lc = c;
  }
  rope_insert_run(rlc->rope, p, lc, rl, 0);
  fprintf(stderr, "[M::%s] Built rope in %.3f sec\n", __func__, realtime() - t);
  free(sa);
}

/**********************/
/** *** QUERYING *** **/
/**********************/
int64_t rlc_lf1(const rlcsa_t *rlc, int64_t x, uint8_t c) {
  int64_t cx[6] = {0, 0, 0, 0, 0, 0};
  rope_rank1a(rlc->rope, x + 1, cx);
  return rlc->C[c] + cx[c] - 1;
}

int rlc_extend(const rlcsa_t *rlc, const qint_t *ik, qint_t ok[6],
               int is_back) {
  int64_t tk[6], tl[6];
  int i;
  rope_rank2a(rlc->rope, ik->x[!is_back], ik->x[!is_back] + ik->x[2], tk, tl);
  for (i = 0; i < 6; ++i) {
    ok[i].x[!is_back] = rlc->C[i] + tk[i];
    ok[i].x[2] = (tl[i] -= tk[i]);
  }
  ok[0].x[is_back] = ik->x[is_back];
  ok[4].x[is_back] = ok[0].x[is_back] + tl[0];
  ok[3].x[is_back] = ok[4].x[is_back] + tl[4];
  ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
  ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
  ok[5].x[is_back] = ok[1].x[is_back] + tl[1];
  return 0;
}

/*********************/
/** *** MERGING *** **/
/*********************/
void report_positions(const rlcsa_t *rlc, const uint8_t *seq, uint32_t n,
                      int64_t *positions) {
  int64_t current = rlc->cnts[0] - 1;
  positions[n] = current; // immediately after current
  for (int64_t i = n - 1; i >= 0; --i) {
    uint8_t c = seq[i];
    current = rlc_lf1(rlc, current, c);
    positions[i] = current; // immediately after current
  }
}

int64_t *get_marks(rlcsa_t *rlc, const uint8_t *seq, uint32_t n,
                   uint_kv separators, int nt) {
  uint32_t i, begin, sequences = 0;
  for (i = 0; i < n; ++i) {
    if (seq[i] == 0) {
      kv_push(uint32_t, separators, i);
      ++sequences;
    }
  }

  int64_t *marks = (int64_t *)calloc(
      n, sizeof(int64_t)); // FIXME: warning "exceeds maximum object size"
  uint chunk = max((uint)1, sequences / (8 * nt));
#pragma omp parallel for schedule(dynamic, chunk)
  for (i = 0; i < sequences; i++) {
    // begin = kv_A(end_markers, i);
    begin = i > 0 ? kv_A(separators, i - 1) + 1 : 0;
    report_positions(rlc, seq + begin, kv_A(separators, i) - begin,
                     marks + begin);
  }
  return marks;
}

void merge_ropes(rope_t *rope1, rope_t *rope2, int64_t *marks) {
  rpitr_t *it1 = calloc(1, sizeof(rpitr_t));
  rope_itr_first(rope1, it1);

  // debug
  // uint8_t *b1;
  // while ((b1 = (uint8_t *) rope_itr_next_block(it1)) != 0)
  //     rle_print(b1, 1);

  rpitr_t *it2 = calloc(1, sizeof(rpitr_t));
  rope_itr_first(rope2, it2);
  uint8_t *b2, *b2_i, *b2_end;
  uint8_t c2;
  int64_t l2, m = 0;
  while ((b2 = (uint8_t *)rope_itr_next_block(it2)) != 0) {
    // rle_print(b2, 0);
    b2_i = b2 + 2;
    b2_end = b2 + 2 + *rle_nptr(b2);
    while (b2_i < b2_end) {
      rle_dec1(b2_i, c2, l2);
      while (l2 > 0) {
        // TODO: insert runs instead of single symbols
        rope_insert_run(rope1, marks[m], c2, 1, 0);
        ++m;
        --l2;
      }
    }
  }

  // debug
  // rope_itr_first(rope1, it1);
  // while ((b1 = (uint8_t *) rope_itr_next_block(it1)) != 0)
  //     rle_print(b1, 1);

  free(it1);
  free(it2);
}

void rlc_merge(rlcsa_t *rlc1, rlcsa_t *rlc2, const uint8_t *seq, int nt) {
  double t = realtime();
  uint32_t i;
  uint8_t c;

  // Get $ positions (separators) and "marked positions" (marks)
  uint_kv separators;
  kv_init(separators);
  int64_t *marks = get_marks(rlc1, seq, rlc2->l, separators, nt);
  fprintf(stderr, "[M::%s] Computed marks in %.3f sec..\n", __func__,
          realtime() - t);
  t = realtime();

  // radix_sort(marks, 0, rlc2->l, 24);
  ks_mergesort(int64_t, rlc2->l, marks, 0);
  fprintf(stderr, "[M::%s] Sorted marks in %.3f sec..\n", __func__,
          realtime() - t);

#pragma omp parallel for num_threads(nt) schedule(static)
  for (i = 0; i < rlc2->l; ++i)
    marks[i] += i + 1;

  // Build C
  for (c = 0; c < 6; ++c)
    rlc1->cnts[c] += rlc2->cnts[c];
  rlc1->C[0] = 0;
  for (c = 1; c < 6; ++c)
    rlc1->C[c] = rlc1->C[c - 1] + rlc1->cnts[c - 1];

  rlc1->l += rlc2->l;

  // Merge ropes using "marked positions"
  t = realtime();
  merge_ropes(rlc1->rope, rlc2->rope, marks);
  fprintf(stderr, "[M::%s] Merged ropes in %.3f sec..\n", __func__,
          realtime() - t);

  free(marks);
  free(rlc2);
  kv_destroy(separators);
}
