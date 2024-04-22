#include "rlcsa.hpp"

bool key_comparison(const std::pair<uint32_t, uint32_t> &a, const std::pair<uint32_t, uint32_t> &b) 
{ 
    // Custom comparison logic 
    return a.second < b.second; // it sorts in ascending order 
}

rlcsa_t *rlc_init() {
    rlcsa_t *rlc;
    rlc = (rlcsa_t *) calloc(1, sizeof(rlcsa_t));
    rlc->cnts = (int64_t *) calloc(6, sizeof(int64_t));
    rlc->C = (int64_t *) calloc(6, sizeof(int64_t));
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
  // rlc_merge(rlc, rlc2, sequence, nt);
}

/**************************/
/** *** CONSTRUCTION *** **/
/**************************/
uint32_t setRanks(skew_pair *pairs, uint32_t *keys, std::vector<ss_range> &unsorted, uint nt, uint chunk) {
    std::vector<ss_range> buffers[nt];
    uint subtotals[nt];
    for(uint32_t i = 0; i < nt; i++)
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
            buffers[thread].push_back(std::make_pair(prev - pairs, (curr - 1) - pairs));
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

sa_t *prefixDoubling(sa_t *pairs, uint32_t *keys, std::vector<ss_range> &unsorted, uint32_t n, uint32_t total, int h, int nt) {
    // Double prefix length until sorted
    while (total > 0) {
        uint chunk = std::max((int64_t)1, (int64_t)unsorted.size() / (nt * nt));
#pragma omp parallel for num_threads(nt) schedule(dynamic, chunk)
        for (uint32_t i = 0; i < unsorted.size(); ++i) {
            // Set sort keys for the current range.
            for (sa_t *curr = pairs + unsorted.at(i).first; curr <= pairs + unsorted.at(i).second; ++curr)
                curr->second = keys[curr->first + h];
            std::sort(pairs + unsorted[i].first, pairs + unsorted[i].second + 1, key_comparison);
        }
        total = setRanks(pairs, keys, unsorted, nt, chunk);
        h *= 2;
        // fprintf(stderr, "Sorted with %d, unsorted total = %ld (%ld ranges)\n", h, total, kv_size(*unsorted));
    }

#pragma omp parallel for num_threads(nt) schedule(static)
    for (uint32_t i = 0; i < n; i++)
        pairs[i].second = keys[i];
    return pairs;
}

sa_t *simpleSuffixSort(const uint8_t *sequence, uint32_t n, uint32_t nsep, int nt) {
    if (sequence == 0 || n == 0)
       // TODO: fail more gracefully
       exit(1);

//   double t = realtime();
    std::vector<ss_range> unsorted;
    skew_pair *pairs = (skew_pair *)calloc(n + 1, sizeof(skew_pair));
    uint32_t *keys = (uint32_t *)calloc(n + 1, sizeof(uint32_t)); // In text order.

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
    // fprintf(stderr, "[M::%s] Sorted pairs in %.3f sec\n", __func__, realtime() - t);
    // t = realtime();
    unsorted.push_back(std::make_pair(0, n-1));
    uint32_t total = setRanks(pairs, keys, unsorted, nt, 1);
    // fprintf(stderr, "[M::%s] Set ranks in %.3f sec\n", __func__, realtime() - t);

    // TODO: doubling vs tripling?

    // t = realtime();
    sa_t *sa = prefixDoubling(pairs, keys, unsorted, n, total, h, nt);
    // fprintf(stderr, "[M::%s] Prefix doubling in %.3f sec\n", __func__, realtime() - t);
    // free(keys);
    // kv_destroy(unsorted);
    return sa;
}

void rlc_build(rlcsa_t *rlc, const uint8_t *sequence, uint32_t n, int nt) {
    // double t;
    // uint32_t i, rl, p;
    // uint8_t c, lc;
    rlc->l = (int64_t)n;

    // Build C
    for (uint32_t i = 0; i < n; ++i)
        ++rlc->cnts[sequence[i]];
    rlc->C[0] = 0;
    for (uint8_t c = 1; c < 6; ++c) // FIXME: hardcoded
        rlc->C[c] = rlc->C[c - 1] + rlc->cnts[c - 1];

    // Build SA
    // t = realtime();
    sa_t *sa = simpleSuffixSort(sequence, n, (uint32_t)rlc->cnts[0], nt);
    // fprintf(stderr, "[M::%s] Built SA in %.3f sec\n", __func__, realtime() - t);

    // Print BWT
    for (uint32_t i = 0; i < n; ++i) {
        uint8_t c = sequence[sa[i].first == 0 ? n : sa[i].first - 1];
        printf("%c", "$ACGTN"[c]);
    }
    printf("\n");

    // Build RLCSA
#pragma omp parallel for schedule(dynamic, 1)
    for(uint8_t c = 0; c < 6; ++c) {
        if (rlc->cnts[c] == 0) {
            // this->array[c] = 0;
            continue;
        }
        std::pair<uint32_t, uint32_t> *curr = sa + rlc->C[c];
        std::pair<uint32_t, uint32_t> *limit = curr + rlc->cnts[c];
        
        // PsiVector::Encoder encoder(block_size);
        std::pair<uint32_t, uint32_t> run((*curr).first, 1);
        ++curr;

        for(; curr < limit; ++curr) {
            if((*curr).first == run.first + run.second) {
                run.second++;
            } else {
                // encoder.addRun(run.first, run.second);
                run = std::make_pair((*curr).first, 1);
            }
        }
        // encoder.addRun(run.first, run.second);
        // encoder.flush();
        // this->array[c] = new PsiVector(encoder, rlc->l);
    }
    free(sa);
}

// /*********************/
// /** *** MERGING *** **/
// /*********************/
// void report_positions(const rlcsa_t *rlc, const uint8_t *seq, uint32_t n,
//                       int64_t *positions) {
//   int64_t current = rlc->cnts[0] - 1;
//   positions[n] = current; // immediately after current
//   for (int64_t i = n - 1; i >= 0; --i) {
//     uint8_t c = seq[i];
//     current = rlc_lf1(rlc, current, c);
//     positions[i] = current; // immediately after current
//   }
// }

// int64_t *get_marks(rlcsa_t *rlc, const uint8_t *seq, uint32_t n,
//                    uint_kv separators, int nt) {
//   uint32_t i, begin, sequences = 0;
//   for (i = 0; i < n; ++i) {
//     if (seq[i] == 0) {
//       kv_push(uint32_t, separators, i);
//       ++sequences;
//     }
//   }

//   int64_t *marks = (int64_t *)calloc(
//       n, sizeof(int64_t)); // FIXME: warning "exceeds maximum object size"
//   uint chunk = max((uint)1, sequences / (8 * nt));
// #pragma omp parallel for schedule(dynamic, chunk)
//   for (i = 0; i < sequences; i++) {
//     // begin = kv_A(end_markers, i);
//     begin = i > 0 ? kv_A(separators, i - 1) + 1 : 0;
//     report_positions(rlc, seq + begin, kv_A(separators, i) - begin,
//                      marks + begin);
//   }
//   return marks;
// }

// rope_t *merge_ropes(rope_t *rope1, rope_t *rope2, int64_t *marks, uint32_t nm) {
//   rope_t *rope = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN);

//   rpitr_t *it1 = calloc(1, sizeof(rpitr_t));
//   rpitr_t *it2 = calloc(1, sizeof(rpitr_t));
//   rope_itr_first(rope1, it1);
//   rope_itr_first(rope2, it2);

//   uint8_t *b1, *b1_i, *b1_end;
//   uint8_t *b2, *b2_i, *b2_end;
//   uint8_t c1, c2, cc; // decoded char from ropes + current char we are inserting
//   int64_t l1, l2,
//       cl;         // decoded run lenghts + current length of run we are building
//   int64_t m = 0;  // current mark
//   int64_t p1 = 0; // position on run from rope1
//   int64_t p2 = 0; // position on run from rope2
//   int64_t pn_i = 0; // inserting position on new rope
//   int64_t pn = 0;   // current position on new rope

//   b1 = (uint8_t *)rope_itr_next_block(it1);
//   b1_i = b1 + 2;
//   b1_end = b1 + 2 + *rle_nptr(b1);
//   rle_dec1(b1_i, c1, l1);

//   b2 = (uint8_t *)rope_itr_next_block(it2);
//   b2_i = b2 + 2;
//   b2_end = b2 + 2 + *rle_nptr(b2);
//   rle_dec1(b2_i, c2, l2);

//   if (marks[0] == pn) {
//     cc = c2;
//     cl = 1;
//     ++p2;
//     ++m;
//     ++pn;
//   } else {
//     cc = c1;
//     cl = 1;
//     ++p1;
//     ++pn;
//     if (l1 == p1) {
//       // we reinitialize b1 here if needed since in the next while we iterate
//       // over b1, so we need to have it set to 0 if we consumed rope1
//       if (b1_i == b1_end) {
//         b1 = (uint8_t *)rope_itr_next_block(it1);
//         if (b1 != 0) {
//           b1_i = b1 + 2;
//           b1_end = b1 + 2 + *rle_nptr(b1);
//         }
//       }
//       if (b1 != 0)
//         rle_dec1(b1_i, c1, l1);
//       p1 = 0;
//     }
//   }

//   while (b1 != 0) {
//     if (m < nm && pn == marks[m]) {
//       // we first check if we have consumed current run from rope2, if yes, we
//       // move to next run
//       if (l2 == p2) {
//         if (b2_i == b2_end) {
//           b2 = (uint8_t *)rope_itr_next_block(it2);
//           // we don't have to do any check here since if we are here, we are
//           // sure that we have at least another symbol from rope2
//           b2_i = b2 + 2;
//           b2_end = b2 + 2 + *rle_nptr(b2);
//         }
//         if (b2 != 0)
//           rle_dec1(b2_i, c2, l2);
//         p2 = 0;
//       }
//       if (cc == c2) {
//         cl += 1;
//       } else {
//         // fprintf(stderr, "Run: %ld %c\n", cl, "$ACGTN"[cc]);
//         rope_insert_run(rope, pn_i, cc, cl, NULL);
//         pn_i += cl;
//         cl = 1;
//         cc = c2;
//       }
//       ++p2;
//       ++m;
//       /* if (m == nm) */
//       /*   fprintf(stderr, "Consumed rope2, continuing rope1\n"); */
//       ++pn;
//     } else {
//       if (l1 == p1) {
//         if (b1_i == b1_end) {
//           b1 = (uint8_t *)rope_itr_next_block(it1);
//           if (b1 == 0) {
//             /* if (m == nm) */
//             /*   fprintf(stderr, "Consumed rope1. Ending\n"); */
//             /* else */
//             /*   fprintf(stderr, "Consumed rope1, moving to rope2\n"); */
//             break;
//           }
//           b1_i = b1 + 2;
//           b1_end = b1 + 2 + *rle_nptr(b1);
//         }
//         if (b1 != 0)
//           rle_dec1(b1_i, c1, l1);
//         p1 = 0;
//       }
//       if (cc == c1) {
//         cl += 1;
//       } else {
//         // fprintf(stderr, "Run: %ld %c\n", cl, "$ACGTN"[cc]);
//         rope_insert_run(rope, pn_i, cc, cl, NULL);
//         pn_i += cl;
//         cl = 1;
//         cc = c1;
//       }
//       ++p1;
//       ++pn;
//     }
//   }

//   // consume current run from rope2
//   if (p2 < l2) {
//     assert(m < nm);
//     if (cc == c2) {
//       // fprintf(stderr, "1 Run: %ld %c\n", cl + l2 - p2, "$ACGTN"[c2]);
//       rope_insert_run(rope, pn_i, c2, cl + l2 - p2, NULL);
//       pn_i += cl + l2 - p2;
//     } else {
//       // fprintf(stderr, "2 Run: %ld %c\n", cl, "$ACGTN"[cc]);
//       rope_insert_run(rope, pn_i, cc, cl, NULL);
//       pn_i += cl;
//       // fprintf(stderr, "3 Run: %ld %c\n", l2 - p2, "$ACGTN"[c2]);
//       rope_insert_run(rope, pn_i, c2, l2 - p2, NULL);
//       pn_i += l2 - p2;
//     }
//   } else {
//     assert(m < nm);
//     // fprintf(stderr, "4 Run: %ld %c\n", cl, "$ACGTN"[cc]);
//     rope_insert_run(rope, pn_i, cc, cl, NULL);
//     pn_i += cl;
//   }

//   // consume current block from rope2
//   while (b2 != 0 && b2_i < b2_end) {
//     assert(m < nm);
//     rle_dec1(b2_i, c2, l2);
//     // fprintf(stderr, "5 Run: %ld %c\n", l2, "$ACGTN"[c2]);
//     rope_insert_run(rope, pn_i, c2, l2, NULL);
//     pn_i += l2;
//   }
//   // consume rest of rope2
//   while ((b2 = (uint8_t *)rope_itr_next_block(it2)) != 0) {
//     assert(m < nm);
//     b2_i = b2 + 2;
//     b2_end = b2 + 2 + *rle_nptr(b2);
//     while (b2_i < b2_end) {
//       rle_dec1(b2_i, c2, l2);
//       // fprintf(stderr, "6 Run: %ld %c\n", l2, "$ACGTN"[c2]);
//       rope_insert_run(rope, pn_i, c2, l2, NULL);
//       pn_i += l2;
//     }
//   }
//   /* fprintf(stderr, "Length: %ld\n", pn_i); */
//   free(it1);
//   free(it2);

//   return rope;
// }

// void rlc_merge(rlcsa_t *rlc1, rlcsa_t *rlc2, const uint8_t *seq, int nt) {
//   double t = realtime();
//   uint32_t i;
//   uint8_t c;

//   // Get $ positions (separators) and "marked positions" (marks)
//   uint_kv separators;
//   kv_init(separators);
//   int64_t *marks = get_marks(rlc1, seq, rlc2->l, separators, nt);
//   fprintf(stderr, "[M::%s] Computed marks in %.3f sec..\n", __func__,
//           realtime() - t);
//   t = realtime();

//   // radix_sort(marks, 0, rlc2->l, 24);
//   ks_mergesort(int64_t, rlc2->l, marks, 0);
//   fprintf(stderr, "[M::%s] Sorted marks in %.3f sec..\n", __func__,
//           realtime() - t);

// #pragma omp parallel for num_threads(nt) schedule(static)
//   for (i = 0; i < rlc2->l; ++i)
//     marks[i] += i + 1;

//   // Build C
//   for (c = 0; c < 6; ++c)
//     rlc1->cnts[c] += rlc2->cnts[c];
//   rlc1->C[0] = 0;
//   for (c = 1; c < 6; ++c)
//     rlc1->C[c] = rlc1->C[c - 1] + rlc1->cnts[c - 1];

//   rlc1->l += rlc2->l;

//   // Merge ropes using "marked positions"
//   t = realtime();
//   rope_t *new_rope = merge_ropes(rlc1->rope, rlc2->rope, marks, rlc2->l);
//   rope_destroy(rlc1->rope);
//   rlc1->rope = new_rope;

//   fprintf(stderr, "[M::%s] Merged ropes in %.3f sec..\n", __func__,
//           realtime() - t);

//   free(marks);
//   free(rlc2);
//   kv_destroy(separators);
// }
