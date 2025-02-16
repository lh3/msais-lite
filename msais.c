/*
 * Copyright (c) 2008  Yuta Mori <yuta.256@gmail.com>
 *               2011- Attractive Chaos <attractor@live.co.uk>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/* This library constructs the generalized suffix array for a string set. It is
 * modified from an early version of sais-lite written by Yuta Mori in 2008. */

#include <stdlib.h>
#include "msais.h"

#if defined(_KSA64) || defined(MSAIS64)
typedef int64_t saint_t;
#define SAIS_MAIN ksa_sa64
#else
typedef int32_t saint_t;
#define SAIS_MAIN ksa_sa32
#endif

// T is of type "const uint8_t*"
#define chr0(i) (cs == sizeof(saint_t) ? ((const saint_t *)T)[i] : T[i])

// count the occurrences of each symbol
static void getCounts(const uint8_t *T, saint_t *C, saint_t n, saint_t k, int cs)
{
	saint_t i;
	for (i = 0; i < k; ++i) C[i] = 0;
	for (i = 0; i < n; ++i) ++C[chr0(i)];
}

// find the start or the end of each bucket, depending on _is_end_
static inline void getBuckets(const saint_t *C, saint_t *B, saint_t k, int is_end)
{
	saint_t i, sum = 0;
	if (is_end) for (i = 0; i < k; ++i) sum += C[i], B[i] = sum;
	else for (i = 0; i < k; ++i) sum += C[i], B[i] = sum - C[i]; // NB: don't change because C and B may point to the same address
}

// place LMS at the end of each bucket
static void placeLMS(const uint8_t *T, saint_t *SA, saint_t *C, saint_t *B, saint_t n, saint_t k, int cs)
{
	saint_t i, c, c0, c1;
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	// find ends of buckets
	for (i = 0; i < n; ++i) SA[i] = 0;
	for (i = n - 2, c = 1, c1 = chr0(n - 1); 0 <= i; --i, c1 = c0) {
		if ((c0 = chr0(i)) < c1 + c) c = 1; // c1 = chr(i+1); c==1 if in an S run
		else if (c) SA[--B[c1]] = i + 1, c = 0;
	}
}

// induce LML from LMS; LMS needs to be correctly placed
static void induceLML(const uint8_t *T, saint_t *SA, saint_t *C, saint_t *B, saint_t n, saint_t k, int cs)
{
	saint_t *b, i, j, c0, c1;
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 0);	// find starts of buckets
	for (i = 0, b = SA, c1 = 0; i < n; ++i) {
		j = SA[i], SA[i] = ~j;
		if (j > 0) { // L or LMS
			--j;
			if ((c0 = chr0(j)) != c1) // then change a bucket
				B[c1] = b - SA, b = SA + B[c1 = c0];
			*b++ = j > 0 && chr0(j - 1) < c1? ~j : j; // true if j is LML, which is <0 now but will be flipped later
		}
	} // at the end of the loop, only LML are positive in SA[]
}

// induce all entries if _all_ is true; otherwise, induce LMS only and set the rest to 0 or negative
static void induceS(const uint8_t *T, saint_t *SA, saint_t *C, saint_t *B, saint_t n, saint_t k, int cs, int all)
{
	saint_t *b, i, j, c0, c1;
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	// find ends of buckets
	if (!all) // set non-LML entries to 0 if we want to induce LMS only
		for (i = B[0]; i < n; ++i)
			SA[i] = SA[i] >= 0? SA[i] : 0;
	for (i = n - 1, b = SA + B[c1 = 0]; i >= 0; --i) {
		j = SA[i];
		if (!all || j <= 0) SA[i] = ~j;
		if (j > 0) {
			--j;
			if ((c0 = chr0(j)) != c1)
				B[c1] = b - SA, b = SA + B[c1 = c0];
			if (c0 > 0) // don't touch the 0 bucket
				*--b = j == 0 || chr0(j - 1) > c1? ~j : j; // true if j is LMS
		}
	}
}

static saint_t nameLMS(const uint8_t *T, saint_t *SA, saint_t n, saint_t m, saint_t k, int cs)
{
	saint_t i, j, name, q, qsig, c, c0, c1, x = k, nb = 0, max_len, n0, shift = sizeof(saint_t) * 8 - 1;
	// calculate max length
	while (x >>= 1) ++nb; // nb = floor(log2(k))
	if ((k & (k - 1)) != 0) ++nb; // nb = ceil(log2(k))
	max_len = shift / nb;
	// store the length or the actual content of all LMS-substrings
	for (i = m; i < n; ++i) SA[i] = 0;
	for (i = n - 2, j = n, c = 1, c1 = chr0(n - 1), x = 0, n0 = 0; i >= 0; --i) {
		c0 = chr0(i);
		if (c0 < c1 + c) c = 1;
		else if (c) {
			saint_t len = j - i - 1;
			SA[m + ((i + 1) >> 1)] = len > max_len || n0 > 0? len : (saint_t)1<<shift | x;
			x = 0, n0 = 0, j = i + 1, c = 0;
		}
		x = x << nb | c0; // substring since the last LMS
		n0 += (c0 == 0); // the number of sentinels since the last LMS
		c1 = c0;
	}
	// find the lexicographic names of all LMS-substrings
	for (i = 0, name = 0, q = n, qsig = 0; i < m; ++i) {
		saint_t p = SA[i], psig = SA[m + (p >> 1)], diff = 1;
		if (psig == qsig) {
			if (psig >> shift == 0) { // then psig is the length of the LMS-substring
				saint_t len = psig;
				for (j = 0; j < len; j++) {
					c0 = chr0(p + j), c1 = chr0(q + j);
					if (c0 != c1 || c0 == 0) break;
				}
				if (j == len) diff = 0;
			} else diff = 0;
		}
		if (diff) ++name, q = p, qsig = psig;
		SA[m + (p >> 1)] = name; // this works because a) m <= n/2 and b) an LMS-substring is at least 2 in length
	}
	return name;
}

/**
 * Recursively construct the suffix array for a string containing multiple
 * sentinels. NULL is taken as the sentinel.
 *
 * @param T   NULL terminated input string (there can be multiple NULLs)
 * @param SA  output suffix array
 * @param fs  working space available in SA (typically 0 when first called)
 * @param n   length of T, including the trailing NULL
 * @param k   size of the alphabet (typically 256 when first called)
 * @param cs  bytes per symbol; typically 1 for the first iteration
 *
 * @return    0 upon success
 */
static int sais_core(const uint8_t *T, saint_t *SA, saint_t fs, saint_t n, saint_t k, int cs)
{
	saint_t *C, *B;
	saint_t  i, j, c, c0, c1, m, max_name;

	// STAGE I: reduce the problem by at least 1/2 sort all the S-substrings
	if (k <= fs) C = SA + n, B = (k <= fs - k) ? C + k : C;
	else {
		if ((C = (saint_t*)malloc(k * (1 + (cs == 1)) * sizeof(saint_t))) == NULL) return -2;
		B = cs == 1? C + k : C;
	}
	placeLMS(T, SA, C, B, n, k, cs);
	induceLML(T, SA, C, B, n, k, cs);
	induceS(T, SA, C, B, n, k, cs, 0);
	if (fs < k) free(C);
	for (i = 0, m = 0; i < n; ++i) // gather all LMS
		if (SA[i] > 0) SA[m++] = SA[i];
	max_name = nameLMS(T, SA, n, m, k, cs);

	// STAGE II: solve the reduced problem; recurse if names are not yet unique
	if (max_name < m) {
		saint_t *RA = SA + n + fs - m - 1; // RA points to the last m+1 elements in SA
		for (i = n - 1, j = m - 1; i >= m; --i)
			if (SA[i] != 0) RA[j--] = SA[i];
		RA[m] = 0; // add a sentinel; in the resulting SA, SA[0]==m
		if (sais_core((uint8_t*)RA, SA, fs + n - m * 2 - 2, m + 1, max_name + 1, sizeof(saint_t)) != 0) return -2;
		for (i = n - 2, j = m - 1, c = 1, c1 = chr0(n - 1); 0 <= i; --i, c1 = c0) {
			if ((c0 = chr0(i)) < c1 + c) c = 1;
			else if (c) RA[j--] = i + 1, c = 0;
		}
		for (i = 0; i < m; ++i) SA[i] = RA[SA[i+1]];
	}

	// STAGE III: induce the result for the original problem
	if (k <= fs) C = SA + n, B = (k <= fs - k) ? C + k : C;
	else {
		if ((C = (saint_t*)malloc(k * (1 + (cs == 1)) * sizeof(saint_t))) == NULL) return -2;
		B = cs == 1? C + k : C;
	}
	// put all LMS characters into their buckets
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	// find ends of buckets
	for (i = m; i < n; ++i) SA[i] = 0;
	for (i = m - 1; i >= 0; --i) {
		j = SA[i], SA[i] = 0;
		SA[--B[chr0(j)]] = j;
	}
	induceLML(T, SA, C, B, n, k, cs);
	induceS(T, SA, C, B, n, k, cs, 1);
	if (fs < k) free(C);
	return 0;
}

/**
 * Construct the suffix array for a NULL terminated string possibly containing
 * multiple sentinels (NULLs).
 *
 * @param T[0..n-1]  NULL terminated input string
 * @param SA[0..n-1] output suffix array
 * @param n          length of the given string, including NULL
 * @param k          size of the alphabet including the sentinel; no more than 256
 * @return           0 upon success
 */
int SAIS_MAIN(const uint8_t *T, saint_t *SA, saint_t n, int k)
{
	if (T == NULL || SA == NULL || n <= 0 || T[n - 1] != '\0') return -1;
	if (k < 0 || k > 256) k = 256;
	return sais_core(T, SA, 0, n, (saint_t)k, 1);
}
