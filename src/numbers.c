#include "numbers.h"

#include <stdio.h>

typedef union 
{
    double fN;
    u64    iN;
} IntFloatToFloatInt;

inline u64 GetFloatNumberBytes(double number)
{
      IntFloatToFloatInt reInterpret = {0};
      reInterpret.fN = number;

      return reInterpret.iN;
}

inline u64 GetFloatNumberMantissa(const u64 numberBytes)
{
      return (numberBytes & 0x000FFFFFFFFFFFFF);
}

inline i32 GetFloatNumberBiasExponent(const u64 numberBytes)
{
      return (i32) ((numberBytes >> 52) & 0x00000000000007FF);
}

inline i32 GetFloatNumberExponent(const u64 numberBytes)
{
      return (i32) ((numberBytes >> 52) & 0x00000000000007FF) - 1023;
}

inline i32 GetFloatNumberSignal(const u64 numberBytes)
{
      return (i32) ((numberBytes >> 63) & 0x0000000000000001);
}

inline double CreateFloatNumber(const u64 mantissa, const u32 biasExp, const u32 signal)
{
      u64 numberBytes = 
      (((u64) (signal & 0x00000001))) << 63 | 
      ((u64) ((biasExp & 0x000007FF))) << 52 | 
      (mantissa & 0x000FFFFFFFFFFFFFLLU);
      
      IntFloatToFloatInt reInterpret = {0};
      reInterpret.iN = numberBytes;

      return reInterpret.fN;
}

void NormalizeFloatNumber(double number, double* pNormalized, i32* pExponent, i32* pSignal)
{
      const u64 numberBytes = GetFloatNumberBytes(number);
      int hasCorrection = 0;

      u64 numberMantissa = GetFloatNumberMantissa(numberBytes);
      i32 numberSignal = GetFloatNumberSignal(numberBytes);
      i32 numberExponent = GetFloatNumberExponent(numberBytes);

      if ((numberExponent % 2) != 0)
      {
            numberExponent -= 1;
            hasCorrection = 1;
      }

      double normalizedNumber = CreateFloatNumber(numberMantissa, 1023, 0);

      *pNormalized = normalizedNumber * (hasCorrection ? 2.0 : 1.0);
      *pExponent = numberExponent;
      *pSignal = numberSignal;
}