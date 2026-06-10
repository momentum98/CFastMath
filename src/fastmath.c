#include "fastmath.h"

#include <immintrin.h>

typedef union
{
      float fN;
      u64   uN;
} U64FloatToFloatU64;

typedef union
{
      float fN;
      u32   uN;
} U32FloatToFloatU32;

typedef union
{
      double fN;
      u64    uN;
} U64DoubleToDoubleU64;

typedef union
{
      double fN;
      u32    uN;
} U32DoubleToDoubleU32;

double SoftwareSqrtFrstLvl(double number, int iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      U64DoubleToDoubleU64 reInterpret = {0};

      reInterpret.fN = number;

      const u64 numberBytes = reInterpret.uN;
            
      const u32 biasExponent = (numberBytes >> 52) & 0x00000000000007FF;
      i32 numberExponent = (i32) biasExponent - 1023;

      const i32 isOdd = numberExponent & 1;
            
      numberExponent -= isOdd;
      
      const u64 normalizedBytes = ((u64) (1023 + isOdd) << 52) | (numberBytes & 0x000FFFFFFFFFFFFFLLU);
      
      reInterpret.uN = normalizedBytes;

      double normalizedNumber = reInterpret.fN;
      
      double try = 0.5 * (1.0 + normalizedNumber);
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = 0.5 * (try + (normalizedNumber / try));
      }

      const i32 finalExponent = (numberExponent >> 1) + 1023;

      reInterpret.fN = try;

      const u64 tryBytes = reInterpret.uN;
      const u64 finalBytes = ((u64) finalExponent << 52) | (tryBytes & 0x000FFFFFFFFFFFFFLLU);

      reInterpret.uN = finalBytes;

      return reInterpret.fN;
}

float SoftwareSqrtFrstLvlF(float number, int iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      U32FloatToFloatU32 reInterpret = {0};

      reInterpret.fN = number;

      const u32 numberBytes = reInterpret.uN;
            
      const u32 biasExponent = (numberBytes >> 23) & 0x00000000000000FF;
      i32 numberExponent = (i32) biasExponent - 127;

      const u32 isOdd = numberExponent & 1;
    
      numberExponent -= isOdd;

      const u32 normalizedBytes = ((u32) (127 + isOdd) << 23) | (numberBytes & 0x00000000007FFFFF);

      reInterpret.uN = normalizedBytes;

      float normalizedNumber = reInterpret.fN;

      float try = 0.5 * (1.0 + normalizedNumber);
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = 0.5 * (try + (normalizedNumber / try));
      }

      const i32 finalExponent = (numberExponent >> 1) + 127;

      reInterpret.fN = try;

      const u32 tryBytes = reInterpret.uN;
      const u32 finalBytes = ((u32) finalExponent << 23) | (tryBytes & 0x00000000007FFFFF);

      reInterpret.uN = finalBytes;

      return reInterpret.fN;
}

double SoftwareSqrtScndLvl(double number, int iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      U64DoubleToDoubleU64 reInterpret = {0};

      reInterpret.fN = number;

      const u64 numberBytes = reInterpret.uN;
            
      const u32 biasExponent = (numberBytes >> 52) & 0x00000000000007FF;
      i32 numberExponent = (i32) biasExponent - 1023;

      const i32 isOdd = numberExponent & 1;
            
      numberExponent -= isOdd;
      
      const u64 normalizedBytes = ((u64) (1023 + isOdd) << 52) | (numberBytes & 0x000FFFFFFFFFFFFFLLU);

      reInterpret.uN = normalizedBytes;

      double normalizedNumber = reInterpret.fN;

      double try = 0.5 * (1.0 + normalizedNumber);
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = try * ((try * try) + (3 * normalizedNumber)) / ((3 * (try * try)) + normalizedNumber);
      }

      const i32 finalExponent = (numberExponent >> 1) + 1023;

      reInterpret.fN = try;

      const u64 tryBytes = reInterpret.uN;
      const u64 finalBytes = ((u64) finalExponent << 52) | (tryBytes & 0x000FFFFFFFFFFFFFLLU);

      reInterpret.uN = finalBytes;

      return reInterpret.fN;
}

float SoftwareSqrtScndLvlF(float number, int iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      U32FloatToFloatU32 reInterpret = {};

      reInterpret.fN = number;

      const u32 numberBytes = reInterpret.uN;
            
      const u32 biasExponent = (numberBytes >> 23) & 0x00000000000000FF;
      i32 numberExponent = (i32) biasExponent - 127;

      const u32 isOdd = numberExponent & 1;
    
      numberExponent -= isOdd;

      const u32 normalizedBytes = ((u32) (127 + isOdd) << 23) | (numberBytes & 0x00000000007FFFFF);

      reInterpret.uN = normalizedBytes;

      float normalizedNumber = reInterpret.fN;

      float try = 0.5 * (1.0 + normalizedNumber);
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = try * ((try * try) + (3 * normalizedNumber)) / ((3 * (try * try)) + normalizedNumber);
      }

      const i32 finalExponent = (numberExponent >> 1) + 127;

      reInterpret.fN = try;

      const u32 tryBytes = reInterpret.uN;
      const u32 finalBytes = ((u32) finalExponent << 23) | (tryBytes & 0x00000000007FFFFF);

      reInterpret.uN = finalBytes;

      return reInterpret.fN;
}

double HardwareSqrt(double number)
{
      __m128d reg = _mm_set_sd(number);
      __m128d value = _mm_sqrt_sd(reg, reg);

      return _mm_cvtsd_f64(value);
}

float HardwareSqrtF(float number)
{
      __m128d reg = _mm_set_ss(number);
      __m128d value = _mm_sqrt_ss(reg);

      return _mm_cvtss_f32(value);
}

float HardwareFastInverseSqrtF(float number)
{
      __m128d reg = _mm_set_ss(number);
      __m128d value = _mm_rsqrt_ss(reg);

      return _mm_cvtss_f32(value);
}

float HardwareFastSqrtF(float number)
{
      __m128d reg = _mm_set_ss(number);
      __m128d value = _mm_rsqrt_ss(reg);

      return _mm_cvtss_f32(value) * number;
}