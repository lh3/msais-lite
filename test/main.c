#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "msais.h"
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define Malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define Realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

#define Grow(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = Realloc(type, (ptr), (__m)); \
		} \
	} while (0)

unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

void seq_char2nt6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l; ++i)
		s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
}

void seq_revcomp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

uint32_t SA_checksum(int64_t len, const int32_t *s)
{
	uint32_t h = 2166136261U;
	int64_t i;
	for (i = 0; i < len; ++i)
		h ^= s[i], h *= 16777619;
	return h;
}

uint32_t SA_checksum64(int64_t len, const int64_t *s)
{
	uint32_t h = 2166136261U;
	int64_t i;
	for (i = 0; i < len; ++i)
		h ^= s[i], h *= 16777619;
	return h;
}

double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#if defined(__linux__)
	return r.ru_maxrss * 1024;
#elif defined(__APPLE__)
	return r.ru_maxrss;
#endif
}

double realtime(void)
{
	static double realtime0 = -1.0;
	struct timeval tp;
	double t;
	gettimeofday(&tp, NULL);
	t = tp.tv_sec + tp.tv_usec * 1e-6;
	if (realtime0 < 0.0) realtime0 = t;
	return t - realtime0;
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	kseq_t *seq;
	gzFile fp;
	int64_t l = 0, max = 0;
	int32_t c, use64 = 0, add_rev = 0;
	uint32_t checksum = 0;
	uint8_t *s = 0;
	double t_real, t_cpu;

	while ((c = ketopt(&o, argc, argv, 1, "r6", 0)) >= 0) {
		if (c == 'r') add_rev = 1;
		else if (c == '6') use64 = 1;
	}
	if (argc == o.ind) {
		fprintf(stderr, "Usage: msais-test [options] <input.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -r        include reverse complement sequences\n");
		fprintf(stderr, "  -6        use 64-bit integers\n");
		return 1;
	}

	// read FASTA/Q
	t_real = realtime();
	t_cpu = cputime();
	fp = gzopen(argv[o.ind], "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		Grow(uint8_t, s, l + (seq->seq.l + 2), max); // +2 to leave room for gSACA-K
		seq_char2nt6(seq->seq.l, (uint8_t*)seq->seq.s);
		memcpy(s + l, seq->seq.s, seq->seq.l + 1); // NB: we are copying 0
		l += seq->seq.l + 1;
		if (add_rev) {
			Grow(uint8_t, s, l + (seq->seq.l + 2), max);
			seq_revcomp6(seq->seq.l, (uint8_t*)seq->seq.s);
			memcpy(s + l, seq->seq.s, seq->seq.l + 1);
			l += seq->seq.l + 1;
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	printf("(MM) Read file in %.3f*%.3f sec (Peak RSS: %.3f MB)\n", realtime() - t_real, (cputime() - t_cpu) / (realtime() - t_real), peakrss() / 1024.0 / 1024.0);

	t_real = realtime();
	t_cpu = cputime();
	if (use64) {
		int64_t *SA = Malloc(int64_t, l);
		ksa_sa64(s, SA, l, 6);
		checksum = SA_checksum64(l, SA);
		free(SA); free(s);
	} else { // ksa32
		int32_t *SA = Malloc(int32_t, l);
		ksa_sa32(s, SA, l, 6);
		checksum = SA_checksum(l, SA);
		free(SA); free(s);
	}
	printf("(MM) Generated SA in %.3f*%.3f sec (Peak RSS: %.3f MB; checksum: %x)\n", realtime() - t_real, (cputime() - t_cpu) / (realtime() - t_real), peakrss() / 1024.0 / 1024.0, checksum);
	return 0;
}

