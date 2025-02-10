## Getting Started
```sh
git clone https://github.com/lh3/msais-lite/
cd msais-lite/test
make
./msais-test test.fa
```

## Introduction

msais-lite is a C library for constructing the [generalized suffix array][gsa] of a string set.
It is adapted from an early version of sais-lite written by [Yuta Mori][yuta] in 2008.
The library only has one API:
```c
int ksa_sa32(const uint8_t *T, int32_t *SA, int32_t n, int k);
```
where `T` is the input string in which each 0 is taken as a sentinel,
`n` is its length and `k` is the maximal symbol in `T` plus 1.
Sentinels in `T` are ordered by their positions:
if $`i<j`$ and $`T[i]`$ and $`T[j]`$ are both 0, msais-lite considers $`T[i]<T[j]`$.
msais-lite is faster than [gSACA-K][gsacak] but slower [libsais][libsais].
It is much simpler than both.

## Algorithm

The version of sais-lite this repo is based on differs from
the [original SA-IS implementation][nong-download] by Ge Nong in that:

1. The original SA-IS uses a bit array to keep the L/S-type, while sais-lite computes the type on the fly.
   This reduces peak memory and random memory accesses.
   Although sais-lite in theory takes $`O(n)`$ working space in the worst case,
   on practical inputs, the working space is much smaller than that.

2. The original SA-IS induces S-types from L-types,
   while sais-lite induces from left-most L-types.
   This is faster as other L-types are irrelevant.

To construct generalized suffix array, msais-lite implicitly converts `T[i]=0`
to `T[i]=i-INT32_MAX` on the fly during character comparisons.
This trick correctly orders sentinels and makes sure they are smaller than all other symbols.
Furthermore, learning from [Timoshevskaya and Feng (2014)][sais-opt],
msais-lite only computes the positions of LMS-suffixes when inducing S-types.
It generates packed the array of LMS positions without using `T`.
This reduces random memory access and speeds up construction a little.

## History of msais-lite

[Ko and Aluru (2005)][ko-2005] defined L-suffix and S-suffix.
They could sort L-suffixes if S-suffixes are already ordered, and vice versa.
[Nong et al (2009)][nong-2009] refined the earlier algorithm.
They introduced left-most S-suffix (LMS-suffix) and found
that we can sort L-suffixes if we know the order of LMS-suffixes, a subset of S-suffixes.
Such sorting is named as *induced sort* in the paper.
The authors called their algorithm SA-IS and provided a succinct implementation in C.
[Yuta Mori][yuta], who developed [libdivsufsort][libdivsufsort],
optimized the original SA-IS implementation in 2008 with improvements shown above
and released it as the *sais-lite* library.
I copied this 2008 version [to bwa][bwa-is].
[Timoshevskaya and Feng (2014)][sais-opt] called it *bwa-is*, but this is Yuta's work, not mine.
Yuta kept improving sais-lite but then the original website is gone.
According to [Web Archive][archive], v2.4.1 is the last version of sais-lite, released on August 7, 2010.
Someone put a copy of this version [in his GitHub repo][filip] and I [forked it][filip-fork] for backup.
sais-lite is an engineering feat.

In 2011, I needed an algorithm to construct generalized suffix array for my [fermi assembler][fermi].
I modified the old sais-lite [in bwa][bwa-is] with source code [here][ksa].
I didn't explain how it works in the [fermi paper][fermi-paper] and I have not touched the source code in 14 years.
During this time, [Louza et al (2017)][gsacak-paper] published gSACA-K.
It is a popular library for constructing generalized suffix arrays albeit slower than my sais-lite adaptation.

Fast forward in 2025, when I was reading a paper using LMS-substrings,
I realized I have completely forgotten how SA-IS works.
I then reread previous papers, came back to my old code, polished it and made it available in this repo.
msais-lite is not the fastest but it is possibly the smallest with good enough performance.

[gsa]: https://en.wikipedia.org/wiki/Generalized_suffix_array
[yuta]: https://github.com/y-256
[gsacak]: https://github.com/felipelouza/gsa-is
[gsacak-paper]: https://www.sciencedirect.com/science/article/pii/S0304397517302621
[libsais]: https://github.com/IlyaGrebnov/libsais
[nong-download]: https://code.google.com/archive/p/ge-nong/downloads
[sais-opt]: https://ieeexplore.ieee.org/document/6863917
[nong-2009]: https://ieeexplore.ieee.org/document/4976463
[ko-2005]: https://www.sciencedirect.com/science/article/pii/S1570866704000498
[libdivsufsort]: https://github.com/y-256/libdivsufsort
[bwa-is]: https://github.com/lh3/bwa/blob/master/is.c
[archive]: https://web.archive.org/web/20151023010453/https://sites.google.com/site/yuta256/
[filip]: https://github.com/fpopic/bioinf/tree/master/install
[filip-fork]: https://github.com/lh3/bioinf-sais-lite/tree/master/install
[fermi]: https://github.com/lh3/fermi
[ksa]: https://github.com/lh3/fermi/blob/master/ksa.c
[fermi-paper]: https://academic.oup.com/bioinformatics/article/28/14/1838/218887
