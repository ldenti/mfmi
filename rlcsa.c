#include "rlcsa.h"

KSORT_INIT(pair, pair_t, pair_lt)

KSORT_INIT_GENERIC(int64_t)

int64_t setRanks(skew_pair *pairs, int64_t *keys, ss_ranges *unsorted, int nt,
                 int chunk) {
  int64_t i, p;
  ss_ranges *buffers = calloc(nt, sizeof(ss_ranges));
  for (i = 0; i < nt; ++i)
    kv_init(buffers[i]);
  uint64_t *subtotals = calloc(nt, sizeof(uint64_t));

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

sa_t *prefixDoubling(sa_t *pairs, int64_t *keys, ss_ranges *unsorted,
                     uint64_t n, int64_t total, int h, int nt) {
  // Double prefix length until sorted
  int64_t i;
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

sa_t *simpleSuffixSort(const uint8_t *sequence, int64_t n, int64_t nsep,
                       int nt) {
  double t = realtime();
  if (sequence == 0 || n == 0)
    // TODO: fail more gracefully
    exit(1);
  ss_ranges unsorted;
  kv_init(unsorted);
  skew_pair *pairs = (skew_pair *)calloc(n + 1, sizeof(skew_pair));
  int64_t *keys = (int64_t *)calloc(n + 1, sizeof(int64_t)); // In text order.

  // Determine pack factor
  int64_t alphabet_size = nsep + 5;          // FIXME: hardcoded
  int64_t limit = INT32_MAX / alphabet_size; // FIXME
  uint h = 1, pack_multiplier = 1;
  while (pack_multiplier * alphabet_size <= limit) {
    ++h;
    pack_multiplier *= alphabet_size;
  }

  // Initialize pairs
  uint zeros = 0, value = 0;
  for (uint i = 0; i < h; ++i) {
    value *= alphabet_size;
    if (sequence[i] == 0) {
      value += zeros;
      zeros++;
    } else {
      value += nsep + sequence[i] - 1;
    }
  }
  for (uint i = 0; i < n - h; i++) {
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
  for (uint i = n - h; i < n; i++) {
    pairs[i].a = i;
    pairs[i].b = value;
    value = (value % pack_multiplier) * alphabet_size;
  }
  // fprintf(stderr, "[M::%s] Built %ld pairs in %.3f sec..\n", __func__, n,
  // realtime() - t); t = realtime();

  // Sort according to first character
  ks_mergesort(pair, n, pairs, 0);
  // fprintf(stderr, "[M::%s] Sorted pairs in %.3f sec..\n", __func__,
  // realtime() - t); t = realtime();
  ss_range a = (ss_range){0, n - 1};
  kv_push(ss_range, unsorted, a);
  int64_t total = setRanks(pairs, keys, &unsorted, nt, 1);
  // fprintf(stderr, "[M::%s] Set ranks in %.3f sec..\n", __func__, realtime() -
  // t);

  // TODO: doubling vs tripling?

  t = realtime();
  sa_t *sa = prefixDoubling(pairs, keys, &unsorted, n, total, h, nt);
  fprintf(stderr, "[M::%s] Prefix doubling in %.3f sec..\n", __func__,
          realtime() - t);
  free(keys);
  kv_destroy(unsorted);
  return sa;
}

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

int rlc_dump(rlcsa_t *rlc, const char *fn) {
  FILE *fp = strcmp(fn, "-") ? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
  if (fp == 0)
    return -1;
  fwrite("RLC\3", 1, 4, fp); // write magic
  fwrite(&rlc->l, 8, 1, fp);
  fwrite(rlc->C, 8, 6, fp);
  fwrite(rlc->cnts, 8, 6, fp);
  rope_dump(rlc->rope, fp);
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

void rlc_build(rlcsa_t *rlc, const uint8_t *sequence, uint32_t n, int nt) {
  double t = realtime();
  uint32_t i;
  int64_t rl, p;
  uint8_t c, lc;
  rlc->l = (int64_t)n;

  // Build C
  for (i = 0; i < n; ++i)
    ++rlc->cnts[sequence[i]];
  rlc->C[0] = 0;
  for (c = 1; c < 6; ++c) // FIXME: hardcoded
    rlc->C[c] = rlc->C[c - 1] + rlc->cnts[c - 1];
  fprintf(stderr, "[M::%s] Built C array in %.3f sec..\n", __func__,
          realtime() - t);

  // Build SA
  t = realtime();
  sa_t *sa = simpleSuffixSort(sequence, n, rlc->cnts[0], nt);
  fprintf(stderr, "[M::%s] Built SA in %.3f sec..\n", __func__, realtime() - t);

  t = realtime();
  lc = sequence[sa[0].a == n ? 0 : sa[0].a - 1];
  rl = 1;
  p = 0;
  for (i = 1; i < n; ++i) {
    c = sequence[sa[i].a == n ? 0 : sa[i].a - 1];
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
  fprintf(stderr, "[M::%s] Built rope in %.3f sec..\n", __func__,
          realtime() - t);
  free(sa);
}

void rlc_insert(rlcsa_t *rlc, const uint8_t *sequence, int64_t n, int nt) {
  if (rlc->l == 0) {
    rlc_build(rlc, sequence, n, nt);
    return;
  }
  rlcsa_t *rlc2 = rlc_init();
  rlc_build(rlc2, sequence, n, nt);
  rlc_merge(rlc, rlc2, sequence, nt);
}

sa_t rlc_init_interval(rlcsa_t *rlc, uint8_t c) {
  return (sa_t){rlc->C[c], rlc->C[c] + rlc->cnts[c] - 1};
}

sa_t rlc_lf(rlcsa_t *rlc, sa_t range, uint8_t c) {
  // FIXME: move these into rlcsa_t?
  int64_t cx[6] = {0, 0, 0, 0, 0, 0};
  int64_t cy[6] = {0, 0, 0, 0, 0, 0};
  rope_rank2a(rlc->rope, range.a, range.b + 1, cx, cy);
  return (sa_t){rlc->C[c] + cx[c], rlc->C[c] + cy[c] - 1};
}

int64_t rlc_lf1(const rlcsa_t *rlc, int64_t x, uint8_t c) {
  int64_t cx[6] = {0, 0, 0, 0, 0, 0};
  rope_rank1a(rlc->rope, x + 1, cx);
  return rlc->C[c] + cx[c] - 1;
}

bisa_t flip(bisa_t range) { return (bisa_t){range.rx, range.x, range.l}; }

bisa_t rlc_init_biinterval(rlcsa_t *rlc, uint8_t c) {
  sa_t x = rlc_init_interval(rlc, c);
  sa_t rx = rlc_init_interval(rlc, (c >= 1 && c <= 4) ? 5 - c : c);
  assert(x.b - x.a == rx.b - rx.a);
  return (bisa_t){x.a, rx.a, x.b - x.a + 1};
}

bisa_t rlc_bilf(rlcsa_t *rlc, bisa_t range, uint8_t c, uint8_t backward) {
  // FIXME: improve this
  if (backward) {
    sa_t *ranges = (sa_t *)calloc(6, sizeof(sa_t));
    int64_t cx[6] = {0, 0, 0, 0, 0, 0};
    int64_t cy[6] = {0, 0, 0, 0, 0, 0};
    for (int c_ = 1; c_ < 5; ++c_) {
      rope_rank2a(rlc->rope, range.x, range.x + range.l, cx, cy);
      ranges[c_] = (sa_t){rlc->C[c_] + cx[c_], cy[c_] - cx[c_]};
      cx[0] = cx[1] = cx[2] = cx[3] = cx[4] = cx[5] = 0;
      cy[0] = cy[1] = cy[2] = cy[3] = cy[4] = cy[5] = 0;
    }
    int64_t *ls = (int64_t *)calloc(6, sizeof(int64_t));
    ls[0] = range.rx + 1;
    ls[4] = ls[0] + ranges[0].b;
    for (int c_ = 3; c_ > 0; --c_)
      ls[c_] = ls[c_ + 1] + ranges[c_ + 1].b;
    bisa_t r = (bisa_t){ranges[c].a, ls[c], ranges[c].b};
    free(ranges);
    free(ls);
    return r;
  } else {
    return flip(rlc_bilf(rlc, flip(range), (c >= 1 && c <= 4) ? 5 - c : c, 1));
  }
}

void report_positions(const rlcsa_t *rlc, const uint8_t *seq, int64_t n,
                      int64_t *positions) {
  uint32_t current = rlc->cnts[0] - 1;
  positions[n] = current; // immediately after current
  for (int64_t i = n - 1; i >= 0; --i) {
    uint8_t c = seq[i];
    current = rlc_lf1(rlc, current, c);
    positions[i] = current; // immediately after current
  }
}

int64_t *get_marks(rlcsa_t *rlc, const uint8_t *seq, uint64_t n,
                   int_kv separators, int nt) {
  uint32_t i, begin, sequences = 0;
  for (i = 0; i < n; ++i) {
    if (seq[i] == 0) {
      kv_push(int64_t, separators, i);
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
  uint i, c;

  // Get $ positions (separators) and "marked positions" (marks)
  int_kv separators;
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
  // fprintf(stderr, "[M::%s] Updated marks in %.3f sec..\n", __func__,
  // realtime() - t);

  // Build C
  for (c = 0; c < 6; ++c)
    rlc1->cnts[c] += rlc2->cnts[c];
  rlc1->C[0] = 0;
  for (c = 1; c < 6; ++c)
    rlc1->C[c] = rlc1->C[c - 1] + rlc1->cnts[c - 1];

  t = realtime();

  // Merge ropes using "marked positions"
  merge_ropes(rlc1->rope, rlc2->rope, marks);

  fprintf(stderr, "[M::%s] Merged ropes in %.3f sec..\n", __func__,
          realtime() - t);

  free(marks);
  free(rlc2);
  kv_destroy(separators);
}
