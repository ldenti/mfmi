#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"

KSEQ_INIT(gzFile, gzread)

static unsigned char seq_nt6_table[128] = {
        0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
        5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

static inline uint kputsn(const char *p, uint l, kstring_t *s) {
    if (s->l + l + 1 >= s->m) {
        char *tmp;
        s->m = s->l + l + 2;
        kroundup32(s->m);
        if ((tmp = (char *) realloc(s->s, s->m)))
            s->s = tmp;
        else
            return EOF;
    }
    memcpy(s->s + s->l, p, l);
    s->l += l;
    s->s[s->l] = 0;
    return l;
}

double cputime() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return (double) r.ru_utime.tv_sec + (double) r.ru_stime.tv_sec +
           1e-6 * (double) (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime() {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return (double) tp.tv_sec + (double) tp.tv_usec * 1e-6;
}

int main(int argc, char *argv[]) {
    (void) argc; // suppress unused parameter warning
    double t_start;
    t_start = realtime();

    char *fa_path = argv[1]; // reference
    char *fq_path = argv[2]; // perfect reads
    int64_t m = (int64_t) (.97 * 10 * 1024 * 1024 * 1024) + 1;

    rlcsa_t *rlc = rlc_init();

    gzFile fp = gzopen(fa_path, "rb");
    kseq_t *ks = kseq_init(fp);
    int l;
    uint8_t *s;
    kstring_t buf = {0, 0, 0};
    int i;
    while ((l = kseq_read(ks)) >= 0) {
        s = (uint8_t *) ks->seq.s;

        // change encoding
        for (i = 0; i < l; ++i)
            s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

        // Add forward to buffer
        kputsn((char *) ks->seq.s, ks->seq.l + 1, &buf);

        // Add reverse to buffer
        for (i = 0; i < (l >> 1); ++i) {
            int tmp = s[l - 1 - i];
            tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
            s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
            s[i] = tmp;
        }
        if (l & 1)
            s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
        kputsn((char *) s, ks->seq.l + 1, &buf);

        if (buf.l >= m) {
            double ct = cputime(), rt = realtime();
            rlc_insert(rlc, (const uint8_t *) buf.s, (int64_t) buf.l);
            fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec, %.3f CPU sec\n",
                    __func__, (long) buf.l, realtime() - rt, cputime() - ct);
            buf.l = 0;
        }
    }
    if (buf.l) {
        double ct = cputime(), rt = realtime();
        rlc_insert(rlc, (const uint8_t *) buf.s, (int64_t) buf.l);
        fprintf(stderr, "[M::%s] inserted %ld symbols in %.3f sec, %.3f CPU sec\n",
                __func__, (long) buf.l, realtime() - rt, cputime() - ct);
        buf.l = 0;
    }
    free(buf.s);
    kseq_destroy(ks);
    gzclose(fp);

    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
            realtime() - t_start, cputime());

    sa_t interval;
    // for (int c = 0; c < 6; ++c) {
    //     interval = rlc_init_interval(rlc, c);
    //     fprintf(stderr, "%d: [%d, %d]\n", c, interval.a, interval.b);
    // }

    int errors = 0;
    fp = gzopen(fq_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
        s = (uint8_t *) ks->seq.s;

        // change encoding
        for (i = 0; i < l; ++i)
            s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
        i = l - 1;
        interval = rlc_init_interval(rlc, s[i]);
        --i;
        for (; i >= 0; --i) {
            interval = rlc_lf(rlc, interval, s[i]);
            if (interval.b < interval.a) {
                errors += 1;
                break;
            }
        }
    }
    printf("Errors: %d\n", errors);

    rlc_destroy(rlc);

    fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
            realtime() - t_start, cputime());

    return 0;
}
