#include "rlcsa.hpp"

bool key_comparison(const std::pair<uint32_t, uint32_t> &a,
                    const std::pair<uint32_t, uint32_t> &b) {
  // Custom comparison logic
  return a.second < b.second; // it sorts in ascending order
}

rlcsa_t *rlc_init() {
  rlcsa_t *rlc;
  rlc = (rlcsa_t *)calloc(1, sizeof(rlcsa_t));
  rlc->cnts = (int64_t *)calloc(6, sizeof(int64_t));
  rlc->C = (int64_t *)calloc(6, sizeof(int64_t));
  for (int c = 0; c < 6; ++c)
    rlc->bits[c] = 0;
  return rlc;
}

void rlc_destroy(rlcsa_t *rlc) {
  free(rlc->cnts);
  free(rlc->C);
  for (int c = 0; c < 6; ++c)
    if (rlc->bits[c])
      delete rlc->bits[c];
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

/**************************/
/** *** CONSTRUCTION *** **/
/**************************/
uint32_t setRanks(skew_pair *pairs, uint32_t *keys,
                  std::vector<ss_range> &unsorted, uint nt, uint chunk) {
  std::vector<ss_range> buffers[nt];
  uint subtotals[nt];
  for (uint32_t i = 0; i < nt; i++)
    subtotals[i] = 0;

#pragma omp parallel for num_threads(nt) schedule(dynamic, chunk)
  for (uint32_t i = 0; i < unsorted.size(); ++i) {
    int thread = omp_get_thread_num();
    skew_pair *prev = pairs + unsorted.at(i).first;
    keys[prev->first] = unsorted.at(i).first;
    skew_pair *limit = pairs + unsorted.at(i).second;
    while (prev < limit) {
      skew_pair *curr = prev + 1;
      if (curr->second != prev->second) {
        ++prev;
        keys[prev->first] = prev - pairs;
        continue;
      }
      keys[curr->first] = prev - pairs;
      for (++curr; curr <= limit && curr->second == prev->second; ++curr)
        keys[curr->first] = prev - pairs;
      buffers[thread].push_back(
          std::make_pair(prev - pairs, (curr - 1) - pairs));
      subtotals[thread] += curr - prev;
      prev = curr;
      if (prev <= limit)
        keys[prev->first] = prev - pairs;
    }
  }
  unsorted.clear();
  uint32_t total = 0;
  for (uint32_t i = 0; i < nt; ++i) {
    unsorted.insert(unsorted.end(), buffers[i].begin(), buffers[i].end());
    total += subtotals[i];
  }
  return total;
}

sa_t *prefixDoubling(sa_t *pairs, uint32_t *keys,
                     std::vector<ss_range> &unsorted, uint32_t n,
                     uint32_t total, int h, int nt) {
  // Double prefix length until sorted
  while (total > 0) {
    uint chunk = std::max((int64_t)1, (int64_t)unsorted.size() / (nt * nt));
#pragma omp parallel for num_threads(nt) schedule(dynamic, chunk)
    for (uint32_t i = 0; i < unsorted.size(); ++i) {
      // Set sort keys for the current range.
      for (sa_t *curr = pairs + unsorted.at(i).first;
           curr <= pairs + unsorted.at(i).second; ++curr)
        curr->second = keys[curr->first + h];
      std::sort(pairs + unsorted[i].first, pairs + unsorted[i].second + 1,
                key_comparison);
    }
    total = setRanks(pairs, keys, unsorted, nt, chunk);
    h *= 2;
    // fprintf(stderr, "Sorted with %d, unsorted total = %ld (%ld ranges)\n", h,
    // total, kv_size(*unsorted));
  }

#pragma omp parallel for num_threads(nt) schedule(static)
  for (uint32_t i = 0; i < n; i++)
    pairs[i].second = keys[i];
  return pairs;
}

sa_t *simpleSuffixSort(const uint8_t *sequence, uint32_t n, uint32_t nsep,
                       int nt) {
  if (sequence == 0 || n == 0)
    // TODO: fail more gracefully
    exit(1);

  // double t = realtime();
  std::vector<ss_range> unsorted;
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
    pairs[i].first = i;
    pairs[i].second = value;
    value = (value % pack_multiplier) * alphabet_size;
    if (sequence[i + h] == 0) {
      value += zeros;
      ++zeros;
    } else {
      value += nsep + sequence[i + h] - 1;
    }
  }
  for (i = n - h; i < n; i++) {
    pairs[i].first = i;
    pairs[i].second = value;
    value = (value % pack_multiplier) * alphabet_size;
  }

  // Sort according to first character
  // t = realtime();
  std::sort(pairs, pairs + n, key_comparison);
  // fprintf(stderr, "[M::%s] Sorted pairs in %.3f sec\n", __func__, realtime()
  // - t); t = realtime();
  unsorted.push_back(std::make_pair(0, n - 1));
  uint32_t total = setRanks(pairs, keys, unsorted, nt, 1);
  // fprintf(stderr, "[M::%s] Set ranks in %.3f sec\n", __func__, realtime() -
  // t);

  // TODO: doubling vs tripling?

  // t = realtime();
  sa_t *sa = prefixDoubling(pairs, keys, unsorted, n, total, h, nt);
  // fprintf(stderr, "[M::%s] Prefix doubling in %.3f sec\n", __func__,
  // realtime() - t); free(keys); kv_destroy(unsorted);
  return sa;
}

void rlc_build(rlcsa_t *rlc, const uint8_t *sequence, uint32_t n, int nt) {
  double t;
  // uint32_t i, rl, p;
  // uint8_t c, lc;
  t = realtime();
  rlc->l = (int64_t)n;

  // Build C
  for (uint32_t i = 0; i < n; ++i)
    ++rlc->cnts[sequence[i]];
  rlc->C[0] = 0;
  for (uint8_t c = 1; c < 6; ++c) // FIXME: hardcoded
    rlc->C[c] = rlc->C[c - 1] + rlc->cnts[c - 1];
  fprintf(stderr, "[M::%s] Built counts in %.3f sec\n", __func__,
          realtime() - t);

  // Build SA
  t = realtime();
  sa_t *sa = simpleSuffixSort(sequence, n, (uint32_t)rlc->cnts[0], nt);
  fprintf(stderr, "[M::%s] Built SA in %.3f sec\n", __func__, realtime() - t);

  // Print BWT
  // for (uint32_t i = 0; i < n; ++i) {
  //   uint8_t c = sequence[sa[i].first == 0 ? n : sa[i].first - 1];
  //   fprintf(stderr, "%c", "$ACGTN"[c]);
  // }
  // fprintf(stderr, "\n");

  // Build Psi.
  t = realtime();
#pragma omp parallel for num_threads(nt) schedule(static)
  for (uint i = 0; i < n; i++)
    sa[i].first = sa[(sa[i].first + 1) % n].second;
  fprintf(stderr, "[M::%s] Built PSI in %.3f sec\n", __func__, realtime() - t);
  // for (int i = 0; i < n; ++i) {
  //   fprintf(stderr, "%d %d\n", sa[i].first, sa[i].second);
  // }

  t = realtime();
  // Build RLCSA
#pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
  for (uint8_t c = 0; c < 6; ++c) {
    if (rlc->cnts[c] == 0) {
      rlc->bits[c] = 0;
      continue;
    }
    std::pair<uint32_t, uint32_t> *curr = sa + rlc->C[c];
    std::pair<uint32_t, uint32_t> *limit = curr + rlc->cnts[c];

    CSA::RLEVector::Encoder encoder(32); // FIXME: hardcoded
    std::pair<uint32_t, uint32_t> run((*curr).first, 1);
    ++curr;

    for (; curr < limit; ++curr) {
      if ((*curr).first == run.first + run.second) {
        run.second++;
      } else {
        encoder.addRun(run.first, run.second);
        run = std::make_pair((*curr).first, 1);
      }
    }
    encoder.addRun(run.first, run.second);
    encoder.flush();
    rlc->bits[c] = new CSA::RLEVector(encoder, rlc->l);
    // CSA::RLEVector::Iterator *iter =
    //     new CSA::RLEVector::Iterator(*rlc->bits[c]);
    // printf("%c: ", "$ACGTN"[c]);
    // for (int i = 0; i < n; ++i)
    //   printf("%d", iter->isSet(i));
    // printf("\n");
    // delete iter;
  }
  fprintf(stderr, "[M::%s] Built RLCSA in %.3f sec\n", __func__,
          realtime() - t);
  free(sa);
}

/*********************/
/** *** MERGING *** **/
/*********************/
void report_positions(const rlcsa_t *rlc, const uint8_t *seq, uint32_t n,
                      int64_t *positions) {
  // { this could be a function
  CSA::RLEVector::Iterator **iters = new CSA::RLEVector::Iterator *[6];
  for (uint8_t i = 0; i < 6; ++i) {
    if (rlc->bits[i] == 0)
      iters[i] = 0;
    else
      iters[i] = new CSA::RLEVector::Iterator(*(rlc->bits[i]));
  }
  // }

  int64_t current = rlc->cnts[0] - 1;
  positions[n] = current; // immediately after current
  for (int64_t i = n - 1; i >= 0; --i) {
    uint8_t c = seq[i];
    current = rlc->C[c] + iters[c]->rank(current) - 1;
    positions[i] = current; // immediately after current
  }
}

int64_t *get_marks(rlcsa_t *rlc, const uint8_t *seq, uint32_t n,
                   std::vector<uint32_t> &separators, int nt) {
  uint32_t i, begin, sequences = 0;
  for (i = 0; i < n; ++i) {
    if (seq[i] == 0) {
      separators.push_back(i);
      ++sequences;
    }
  }

  int64_t *marks = (int64_t *)calloc(
      n, sizeof(int64_t)); // FIXME: warning "exceeds maximum object size"
  uint chunk = std::max((uint)1, sequences / (8 * nt));
#pragma omp parallel for schedule(dynamic, chunk)
  for (i = 0; i < sequences; i++) {
    // begin = kv_A(end_markers, i);
    begin = i > 0 ? separators.at(i - 1) + 1 : 0;
    report_positions(rlc, seq + begin, separators.at(i) - begin, marks + begin);
  }
  return marks;
}

CSA::RLEVector *merge(CSA::RLEVector *first, CSA::RLEVector *second,
                      int64_t *marks, uint32_t n2 /*n*/, uint nn /*size*/,
                      uint block_size) {
  if ((first == 0 && second == 0) || marks == 0)
    return 0;

  CSA::RLEVector::Iterator *first_iter = 0;
  CSA::RLEVector::Iterator *second_iter = 0;

  std::pair<int64_t, int64_t> first_run;
  bool first_finished;
  if (first == 0) {
    first_run = std::make_pair(nn, 0);
    first_finished = true;
  } else {
    first_iter = new CSA::RLEVector::Iterator(*first);
    first_run = first_iter->selectRun(0, nn);
    ++first_run.second;
    first_finished = false;
  }

  uint64_t second_bit;
  if (second == 0) {
    second_bit = n2;
  } else {
    second_iter = new CSA::RLEVector::Iterator(*second);
    second_bit = second_iter->select(0);
  }

  CSA::RLEVector::Encoder encoder(block_size);
  for (uint32_t i = 0; i < n2; ++i) {
    while (!first_finished && first_run.first + i < marks[i]) {
      int64_t bits = std::min(first_run.second, marks[i] - i - first_run.first);
      encoder.addRun(first_run.first + i, bits);
      first_run.first += bits;
      first_run.second -= bits;
      if (first_run.second == 0) {
        if (first_iter->hasNext()) {
          first_run = first_iter->selectNextRun(nn);
          first_run.second++;
        } else {
          first_finished = true;
        }
      }
    }

    if (i == second_bit) {
      // marks[i] is one
      encoder.addBit(marks[i]);
      second_bit = second_iter->selectNext();
    }
  }

  while (!first_finished) {
    encoder.addRun(first_run.first + n2, first_run.second);
    if (first_iter->hasNext()) {
      first_run = first_iter->selectNextRun(nn);
      first_run.second++;
    } else {
      first_finished = true;
    }
  }

  delete first_iter;
  delete second_iter;
  delete first;
  delete second;
  encoder.flush();
  return new CSA::RLEVector(encoder, nn);
}

void rlc_merge(rlcsa_t *rlc1, rlcsa_t *rlc2, const uint8_t *seq, int nt) {
  // double t = realtime();
  uint32_t i;
  uint8_t c;

  // Get $ positions (separators) and "marked positions" (marks)
  std::vector<uint32_t> separators;
  int64_t *marks = get_marks(rlc1, seq, rlc2->l, separators, nt);
  // fprintf(stderr, "[M::%s] Computed marks in %.3f sec..\n", __func__,
  // realtime() - t); t = realtime();

  // fprintf(stderr, "Marks: ");
  // for (i = 0; i < rlc2->l; ++i)
  //   fprintf(stderr, "%ld ", marks[i]);
  // fprintf(stderr, "\n");

  // radix_sort(marks, 0, rlc2->l, 24);
  std::sort(marks, marks + rlc2->l);
  // fprintf(stderr, "[M::%s] Sorted marks in %.3f sec..\n", __func__,
  // realtime() - t);

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

#pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
  for (c = 0; c < 6; ++c) {
    rlc1->bits[c] =
        merge(rlc1->bits[c], rlc2->bits[c], marks, rlc2->l, rlc1->l, 32);
    rlc2->bits[c] = 0;
  }

  // fprintf(stderr, "[M::%s] Merged ropes in %.3f sec..\n", __func__,
  // realtime() - t);

  free(marks);
  rlc_destroy(rlc2);
}

void rlc_print_bits(rlcsa_t *rlc) {
  // { this could be a function
  CSA::RLEVector::Iterator **iters = new CSA::RLEVector::Iterator *[6];
  for (uint8_t i = 0; i < 6; ++i) {
    if (rlc->bits[i] == 0)
      iters[i] = 0;
    else
      iters[i] = new CSA::RLEVector::Iterator(*(rlc->bits[i]));
  }
  // }

  uint8_t c;
  for (c = 0; c < 6; ++c) {
    if (iters[c] == 0)
      continue;
    printf("%c: ", "$ACGTN"[c]);
    for (int64_t i = 0; i < rlc->l; ++i) {
      printf("%d", iters[c]->isSet(i));
    }
    printf("\n");
  }
}

void rlc_dump(rlcsa_t *rlc) {
  // fprintf(stderr, "Total length: %ld\n", rlc->l);

  // { this could be a function
  CSA::RLEVector::Iterator **iters = new CSA::RLEVector::Iterator *[6];
  for (uint8_t i = 0; i < 6; ++i) {
    if (rlc->bits[i] == 0)
      iters[i] = 0;
    else
      iters[i] = new CSA::RLEVector::Iterator(*(rlc->bits[i]));
  }
  // }

  uint8_t c;

  rld_t *e = 0;
  rlditr_t di;
  e = rld_init(6, 3);
  rld_itr_init(e, &di, 0);

  int64_t rl = 1;
  uint8_t rc = 0;
  while (!iters[rc]->isSet(0))
    ++rc;
  // fprintf(stderr, "%c", "$ACGTN"[rc]);
  for (int64_t i = 1; i < rlc->l; ++i) {
    c = 0;
    while (!iters[c]->isSet(i))
      ++c;
    if (c == rc)
      ++rl;
    else {
      rld_enc(e, &di, rl, rc);
      rl = 1;
      rc = c;
    }
    // fprintf(stderr, "%c", "$ACGTN"[c]);
  }
  rld_enc(e, &di, rl, rc);

  // fprintf(stderr, "\n");

  rld_enc_finish(e, &di);
  rld_dump(e, "-");
}
