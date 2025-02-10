#ifndef MSAIS_H
#define MSAIS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Constructing the generalized suffix array for a string set
 *
 * @param T     string with 0 taken as sentinels; T[n-1] MUST BE 0
 * @param SA    suffix array of length n
 * @param n     number of symbols
 * @param k     largest symbol plus 1
 *
 * @return 0 on success and -1 on failure
 */
int ksa_sa32(const uint8_t *T, int32_t *SA, int32_t n, int k);

int ksa_sa64(const uint8_t *T, int64_t *SA, int64_t n, int k);

#ifdef __cplusplus
}
#endif

#endif
