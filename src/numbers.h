#ifndef NUMBERS_H
#define NUMBERS_H

#include <stdint.h>

typedef uint64_t  u64;
typedef uint32_t  u32;
typedef int32_t   i32;
typedef int64_t   i64;

inline u64 GetFloatNumberBytes(double number);
inline u64 GetFloatNumberMantissa(const u64 numberBytes);
inline i32 GetFloatNumberBiasExponent(const u64 numberBytes);
inline i32 GetFloatNumberExponent(const u64 numberBytes);
inline i32 GetFloatNumberSignal(const u64 numberBytes);
inline double CreateFloatNumber(const u64 mantissa, const u32 biasExp, const u32 signal);
void NormalizeFloatNumber(double number, double* pNormalized, i32* pExponent, i32* pSignal);

#endif