#ifndef NUMBERS_H
#define NUMBERS_H

#include <stdint.h>

typedef uint64_t  u64;
typedef uint32_t  u32;
typedef int32_t   i32;
typedef int64_t   i64;

typedef union 
{
    double fN;
    u64    iN;
} IntFloatToFloatInt;

static inline u64 GetFloatNumberBytes(double number)
{
      IntFloatToFloatInt reInterpret = {0};
      reInterpret.fN = number;

      return reInterpret.iN;
}

static inline u64 GetFloatNumberMantissa(const u64 numberBytes)
{
      return (numberBytes & 0x000FFFFFFFFFFFFF);
}

static inline i32 GetFloatNumberBiasExponent(const u64 numberBytes)
{
      return (i32) ((numberBytes >> 52) & 0x00000000000007FF);
}

static inline i32 GetFloatNumberExponent(const u64 numberBytes)
{
      return (i32) ((numberBytes >> 52) & 0x00000000000007FF) - 1023;
}

static inline i32 GetFloatNumberSignal(const u64 numberBytes)
{
      return (i32) ((numberBytes >> 63) & 0x0000000000000001);
}

static inline double CreateFloatNumber(const u64 mantissa, const u32 biasExp, const u32 signal)
{
      u64 numberBytes = 
      (((u64) (signal & 0x00000001))) << 63 | 
      ((u64) ((biasExp & 0x000007FF))) << 52 | 
      (mantissa & 0x000FFFFFFFFFFFFFLLU);
      
      IntFloatToFloatInt reInterpret = {0};
      reInterpret.iN = numberBytes;

      return reInterpret.fN;
}

void NormalizeFloatNumber(double number, double* pNormalized, i32* pExponent, i32* pSignal);

#endif