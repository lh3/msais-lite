#ifndef MSAIS_H
#define MSAIS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int ksa_sa32(const uint8_t *T, int32_t *SA, int32_t n, int k);
int ksa_bwt32(uint8_t *T, int32_t n, int k);
int ksa_core32(const uint8_t *T, int32_t *SA, int32_t fs, int32_t n, int32_t k, int cs);

int ksa_sa64(const uint8_t *T, int64_t *SA, int64_t n, int k);
int ksa_bwt64(uint8_t *T, int64_t n, int k);
int ksa_core64(const uint8_t *T, int64_t *SA, int64_t fs, int64_t n, int64_t k, int cs);

#ifdef __cplusplus
}
#endif

#endif
