#include "rlcsa.h"

KSORT_INIT(pair, pair_t, pair_lt)

KSORT_INIT_GENERIC(int64_t)

int64_t setRanks(skew_pair *pairs, int64_t *keys,
                 ss_ranges *unsorted) {// int64_t n, uint8_t threads, uint32_t chunk) {
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
    for (int64_t i = 0; i < kv_size(*unsorted); ++i) {
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
            a = (ss_range) {prev - pairs, (curr - 1) - pairs};
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

sa_t *prefixDoubling(sa_t *pairs, int64_t *keys, ss_ranges *unsorted, uint64_t n,
                     int64_t total, int h) { //uint threads
    // Double prefix length until sorted
    // TODO: parallelize
    int64_t i;
    while (total > 0) {
        // uint chunk = std::max((size_t)1, unsorted.size() / (threads * threads));
        // #pragma omp parallel for schedule(dynamic, chunk)
        for (i = 0; i < kv_size(*unsorted); ++i) {
            // Set sort keys for the current range.
            for (sa_t *curr = pairs + kv_A(*unsorted, i).a;
                 curr <= pairs + kv_A(*unsorted, i).b; ++curr)
                curr->b = keys[curr->a + h];
            ks_mergesort(pair, kv_A(*unsorted, i).b - kv_A(*unsorted, i).a + 1, pairs + kv_A(*unsorted, i).a, 0);
        }
        total = setRanks(pairs, keys, unsorted);//, threads, chunk);
        h *= 2;
        fprintf(stderr, "Sorted with %d, unsorted total = %ld (%ld ranges)\n", h,
                total, kv_size(*unsorted));
    }
    // TODO: parallelize
    // #pragma omp parallel for schedule(static)
    for (i = 0; i < n; i++)
        pairs[i].b = keys[i];
    return pairs;
}

sa_t *simpleSuffixSort(const uint8_t *sequence, int64_t n) { // uint8_t threads
    if (sequence == 0 || n == 0)
        // TODO: fail more gracefully
        exit(1);
    ss_ranges unsorted;
    kv_init(unsorted);
    skew_pair *pairs = (skew_pair *) calloc(n + 1, sizeof(skew_pair));
    int64_t *keys =
            (int64_t *) calloc(n + 1, sizeof(int64_t));// In text order.
    // Initialize pairs.
    // TODO: parallelize
    // #pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < n; ++i) {
        pairs[i].a = i;
        pairs[i].b = sequence[i];
    }
    // Sort according to first character
    ks_mergesort(pair, n, pairs, 0);
    ss_range a = (ss_range) {0, n - 1};
    kv_push(ss_range, unsorted, a);
    int64_t total = setRanks(pairs, keys, &unsorted);// threads, 1);
    sa_t *sa = prefixDoubling(pairs, keys, &unsorted, n, total, 1); // threads
    // TODO: doubling vs tripling?
    free(keys);
    kv_destroy(unsorted);
    return sa;
}

rlcsa_t *rlc_init() {
    rlcsa_t *rlc;
    rlc = calloc(1, sizeof(rlcsa_t));
    rlc->cnts = calloc(6, sizeof(int64_t));
    rlc->C = calloc(6, sizeof(int64_t));
    for (int c = 0; c < 6; ++c)
        rlc->bits[c] = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN);
    return rlc;
}

void rlc_destroy(rlcsa_t *rlc) {
    free(rlc->cnts);
    free(rlc->C);
    for (int c = 0; c < 6; ++c)
        rope_destroy(rlc->bits[c]);
    free(rlc);
}

int rlc_insert(rlcsa_t *rlc, const uint8_t *sequence, int64_t n) {
    int i, c;
    rlc->l = n;

    // Build C
    for (i = 0; i < n; ++i)
        ++rlc->cnts[sequence[i]];
    rlc->C[0] = 0;
    for (c = 1; c < 6; ++c)// FIXME: hardcoded
        rlc->C[c] = rlc->C[c - 1] + rlc->cnts[c - 1];

    // Build SA
    sa_t *sa = simpleSuffixSort(sequence, n); // threads

    // Build Psi
    // TODO: #pragma omp parallel for schedule(static)
    for (i = 0; i < n; ++i)
        sa[i].a = sa[(sa[i].a + 1) % n].b;

    // Build RLCSA
    // TODO: #pragma omp parallel for schedule(dynamic, 1)
    for (c = 0; c < 6; ++c) {// TODO: hardcoded
        if (rlc->cnts[c] == 0)
            continue;

        sa_t *curr = sa + rlc->C[c];
        sa_t *limit = curr + rlc->cnts[c];
        int64_t last_e = 0;
        int64_t curr_s = curr->a;
        int64_t curr_l = 1;
        ++curr;

        for (; curr < limit; ++curr) {
            if (curr->a == curr_s + curr_l)
                ++curr_l;
            else {
                if ((last_e == 0 && curr_s > 0) || last_e != 0) {
                    rope_insert_run(rlc->bits[c], last_e, 0, curr_s - last_e, 0);
                }
                rope_insert_run(rlc->bits[c], curr_s, 1, curr_l, 0);
                last_e = curr_s + curr_l;
                curr_s = curr->a;
                curr_l = 1;
            }
        }
        if ((last_e == 0 && curr_s > 0) || last_e != 0)
            rope_insert_run(rlc->bits[c], last_e, 0, curr_s - last_e, 0);
        rope_insert_run(rlc->bits[c], curr_s, 1, curr_l, 0);
        curr_s += curr_l;
        if (n - curr_s != 0)
            rope_insert_run(rlc->bits[c], curr_s, 0, n - curr_s, 0);
    }
    free(sa);

    return 0;
}

sa_t rlc_init_interval(rlcsa_t *rlc, uint8_t c) {
    return (sa_t) {rlc->C[c], rlc->C[c] + rlc->cnts[c] - 1};
}

sa_t rlc_lf(rlcsa_t *rlc, sa_t range, uint8_t c) {
    // FIXME: move these into rlcsa_t?
    int64_t cx[6] = {0, 0, 0, 0, 0, 0};
    int64_t cy[6] = {0, 0, 0, 0, 0, 0};
    rope_rank2a(rlc->bits[c], range.a, range.b + 1, cx, cy);
    return (sa_t) {rlc->C[c] + cx[1], rlc->C[c] + cy[1] - 1};
}

int64_t rlc_lf1(const rlcsa_t *rlc, int64_t x, uint8_t c) {
    int64_t cx[6] = {0, 0, 0, 0, 0, 0};
    rope_rank1a(rlc->bits[c], x + 1, cx);
    return rlc->C[c] + cx[1] - 1; // FIXME: + rlc->cnts[0];
}

void report_positions(const rlcsa_t *rlc, const uint8_t *seq, int64_t n, int64_t *positions) {
    uint32_t current = rlc->cnts[0] - 1;
    positions[n] = current; // immediately after current
    for (int64_t i = n - 1; i >= 0; --i) {
        uint8_t c = seq[i];
        // FIXME: assuming to have bitvector for character c
        current = rlc_lf1(rlc, current, c);
        positions[i] = current; // immediately after current
    }
}

int64_t *get_marks(rlcsa_t *rlc, const uint8_t *seq, uint64_t n, int_kv separators) {
    uint32_t i, begin, sequences = 0;
    for (i = 0; i < n; ++i) {
        if (seq[i] == 0) {
            kv_push(int64_t, separators, i);
            ++sequences;
        }
    }

    int64_t *marks = (int64_t *) calloc(n, sizeof(int64_t)); // FIXME: warning "exceeds maximum object size"

    // TODO: parallelize
    // usint chunk = std::max((usint)1, sequences / (8 * this->threads));
    // #pragma omp parallel for schedule(dynamic, chunk)
    for (i = 0; i < sequences; i++) {
        // begin = kv_A(end_markers, i);
        begin = i > 0 ? kv_A(separators, i - 1) + 1 : 0;
        report_positions(rlc, seq + begin, kv_A(separators, i) - begin, marks + begin);
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
    while ((b2 = (uint8_t *) rope_itr_next_block(it2)) != 0) {
        // rle_print(b2, 1);
        b2_i = b2 + 2;
        b2_end = b2 + 2 + *rle_nptr(b2);
        while (b2_i < b2_end) {
            rle_dec1(b2_i, c2, l2);
            // printf("%d[%ld]\n", c2, l2);
            while (l2 > 0) {
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

void rlc_merge(rlcsa_t *rlc1, rlcsa_t *rlc2, const uint8_t *seq) { // uint8_t threads
    uint i, c;

    // Get $ positions (separators) and "marked positions" (marks)
    int_kv separators;
    kv_init(separators);
    int64_t *marks = get_marks(rlc1, seq, rlc2->l, separators);
    // TODO: parallelize
    // radix_sort(marks, 0, rlc2->l, 24);
    ks_mergesort(int64_t, rlc2->l, marks, 0);
    // #pragma omp parallel for schedule(static)
    for (i = 0; i < rlc2->l; ++i)
        marks[i] += i + 1;

    // Build C
    for (c = 0; c < 6; ++c)
        rlc1->cnts[c] += rlc2->cnts[c];
    rlc1->C[0] = 0;
    for (c = 1; c < 6; ++c)
        rlc1->C[c] = rlc1->C[c - 1] + rlc1->cnts[c - 1];

    // Merge bit vectors using "marked positions"
    // TODO: parallelize
    // omp_set_num_threads(threads);
    // #pragma omp parallel for schedule(dynamic, 1)
    for (c = 1; c < 6; ++c)
        merge_ropes(rlc1->bits[c], rlc2->bits[c], marks);
    free(marks);
    free(rlc2);
    kv_destroy(separators);
}
