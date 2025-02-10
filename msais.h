#ifndef MSAIS_H
#define MSAIS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int ksa_sa32(const uint8_t *T, int32_t *SA, int32_t n, int k);
int ksa_sa64(const uint8_t *T, int64_t *SA, int64_t n, int k);

#ifdef __cplusplus
}
#endif

#endif
