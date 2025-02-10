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
#define SAINT_MAX INT64_MAX
#define SAIS_MAIN ksa_sa64
#else
typedef int32_t saint_t;
#define SAINT_MAX INT32_MAX
#define SAIS_MAIN ksa_sa32
#endif

// T is of type "const uint8_t*". If T[i] is a sentinel, chr(i) takes a negative value
#define chr(i) (cs == sizeof(saint_t) ? ((const saint_t *)T)[i] : (T[i]? (saint_t)T[i] : i - SAINT_MAX))
#define chr0(i) (cs == sizeof(saint_t) ? ((const saint_t *)T)[i] : T[i])

/** Count the occurrences of each symbol */
static void getCounts(const uint8_t *T, saint_t *C, saint_t n, saint_t k, int cs)
{
	saint_t i;
	for (i = 0; i < k; ++i) C[i] = 0;
	for (i = 0; i < n; ++i) ++C[chr0(i)];
}

/**
 * Find the start or end of each bucket
 *
 * @param C    occurrences computed by getCounts()
 * @param B    start/end of each bucket (out)
 * @param k    size of alphabet
 * @param end  compute the end of bucket if true; otherwise compute the end
 */
static inline void getBuckets(const saint_t *C, saint_t *B, saint_t k, saint_t end)
{
	saint_t i, sum = 0;
	if (end) for (i = 0; i < k; ++i) sum += C[i], B[i] = sum;
	else for (i = 0; i < k; ++i) sum += C[i], B[i] = sum - C[i]; // NB: don't change because C and B may point to the same address
}

/**
 * Induced sort
 *
 * @param T         string
 * @param SA        suffix array with LMS stored at the ends of buckets
 * @param C         array for counting (no need to compute)
 * @param B         array for bucket offsets (no need to compute)
 * @param n         length of T
 * @param k         size of alphabet
 * @param cs        bytes per symbol; typically 1 for the first iteration
 * @param LMS_only  if false, populate all SA values; otherwise, only LMS positions are positive in SA
 */
static void induceSA(const uint8_t *T, saint_t *SA, saint_t *C, saint_t *B, saint_t n, saint_t k, int cs, int LMS_only)
{
	saint_t *b, i, j;
	saint_t  c0, c1;

	// induce L from LMS (left-to-right)
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

	// induce S from LML (right-to-left)
	if (C == B) getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	// find ends of buckets
	if (LMS_only) // set negative values to 0 except the 0/sentinel bucket
		for (i = B[0]; i < n; ++i)
			if (SA[i] < 0) SA[i] = 0;
	for (i = n - 1, b = SA + B[c1 = 0]; i >= 0; --i) {
		j = SA[i];
		if (LMS_only || j <= 0) SA[i] = ~j;
		if (j > 0) {
			--j;
			if ((c0 = chr0(j)) != c1)
				B[c1] = b - SA, b = SA + B[c1 = c0];
			if (c0 > 0) // don't touch the 0 bucket
				*--b = j == 0 || chr0(j - 1) > c1? ~j : j; // true if j is LMS
		}
	}
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
	saint_t  i, j, c, m, q, qlen, name;
	saint_t  c0, c1;

	// STAGE I: reduce the problem by at least 1/2 sort all the S-substrings
	if (k <= fs) C = SA + n, B = (k <= fs - k) ? C + k : C;
	else {
		if ((C = (saint_t*)malloc(k * (1 + (cs == 1)) * sizeof(saint_t))) == NULL) return -2;
		B = cs == 1? C + k : C;
	}
	getCounts(T, C, n, k, cs);
	getBuckets(C, B, k, 1);	// find ends of buckets
	for (i = 0; i < n; ++i) SA[i] = 0;
	// find LMS and keep their positions in the buckets
	for (i = n - 2, c = 1, c1 = chr0(n - 1); 0 <= i; --i, c1 = c0) {
		if ((c0 = chr0(i)) < c1 + c) c = 1; // c1 = chr(i+1); c==1 if in an S run
		else if (c) SA[--B[c1]] = i + 1, c = 0;
	}
	induceSA(T, SA, C, B, n, k, cs, 1);
	if (fs < k) free(C);
	// pack all the sorted LMS into the first m items of SA; 2*m <= n
	for (i = 0, m = 0; i < n; ++i)
		if (SA[i] > 0) SA[m++] = SA[i];
	for (i = m; i < n; ++i) SA[i] = 0;	// init the name array buffer
	// store the length of all substrings
	for (i = n - 2, j = n, c = 1, c1 = chr0(n - 1); i >= 0; --i, c1 = c0) {
		if ((c0 = chr0(i)) < c1 + c) c = 1; // c1 = chr(i+1)
		else if (c) SA[m + ((i + 1) >> 1)] = j - i - 1, j = i + 1, c = 0;
	}
	// find the lexicographic names of all substrings
	for (i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
		saint_t p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
		if (plen == qlen) {
			for (j = 0; j < plen && chr(p + j) == chr(q + j); j++) {}
			if (j == plen) diff = 0;
		}
		if (diff) ++name, q = p, qlen = plen;
		SA[m + (p >> 1)] = name;
	}

	// STAGE II: solve the reduced problem; recurse if names are not yet unique
	if (name < m) {
		saint_t *RA = SA + n + fs - m - 1;
		for (i = n - 1, j = m - 1; m <= i; --i)
			if (SA[i] != 0) RA[j--] = SA[i];
		RA[m] = 0; // add a sentinel; in the resulting SA, SA[0]==m always stands
		if (sais_core((uint8_t*)RA, SA, fs + n - m * 2 - 2, m + 1, name + 1, sizeof(saint_t)) != 0) return -2;
		for (i = n - 2, j = m - 1, c = 1, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
			if ((c0 = chr(i)) < c1 + c) c = 1;
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
	for (i = m - 1; 0 <= i; --i) {
		j = SA[i], SA[i] = 0;
		SA[--B[chr0(j)]] = j;
	}
	induceSA(T, SA, C, B, n, k, cs, 0);
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
