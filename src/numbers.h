#ifndef NUMBERS_H
#define NUMBERS_H

#include <stdint.h>

typedef uint64_t  u64;
typedef uint32_t  u32;
typedef int32_t   i32;
typedef int64_t   i64;

u64 GetFloatNumberBytes(double number);
u64 GetFloatNumberMantissa(const u64 numberBytes);
i32 GetFloatNumberBiasExponent(const u64 numberBytes);
i32 GetFloatNumberExponent(const u64 numberBytes);
i32 GetFloatNumberSignal(const u64 numberBytes);
double CreateFloatNumber(const u64 mantissa, const u32 biasExp, const u32 signal);
void NormalizeFloatNumber(double number, double* pNormalized, i32* pExponent, i32* pSignal);

#endif