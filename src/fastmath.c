#include "fastmath.h"

#include <immintrin.h>

typedef union
{
      f32 fN;
      u64 uN;
} U64F32ToF32U64;

typedef union
{
      f32 fN;
      u32 uN;
} U32F32ToF32U32;

typedef union
{
      f64 fN;
      u64 uN;
} U64F64ToF64U64;

typedef union
{
      f64 fN;
      u32 uN;
} U32F64ToF64U32;

f64 SSqrt(f64 number, i32 iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      U64F64ToF64U64 reInterpret = {0};

      reInterpret.fN = number;

      const u64 numberBytes = reInterpret.uN;
            
      const u32 biasExponent = (numberBytes >> 52) & 0x00000000000007FF;
      i32 numberExponent = (i32) biasExponent - 1023;

      const i32 isOdd = numberExponent & 1;
            
      numberExponent -= isOdd;
      
      const u64 normalizedBytes = ((u64) (1023 + isOdd) << 52) | (numberBytes & 0x000FFFFFFFFFFFFFLLU);
      
      reInterpret.uN = normalizedBytes;

      f64 normalizedNumber = reInterpret.fN;
      
      f64 try = 0.5 * (1.0 + normalizedNumber);
      
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

f32 SSqrtF(f32 number, i32 iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      U32F32ToF32U32 reInterpret = {0};

      reInterpret.fN = number;

      const u32 numberBytes = reInterpret.uN;
            
      const u32 biasExponent = (numberBytes >> 23) & 0x00000000000000FF;
      i32 numberExponent = (i32) biasExponent - 127;

      const u32 isOdd = numberExponent & 1;
    
      numberExponent -= isOdd;

      const u32 normalizedBytes = ((u32) (127 + isOdd) << 23) | (numberBytes & 0x00000000007FFFFF);

      reInterpret.uN = normalizedBytes;

      f32 normalizedNumber = reInterpret.fN;

      f32 try = 0.5 * (1.0 + normalizedNumber);
      
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

f64 SSqrt2(f64 number, i32 iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      U64F64ToF64U64 reInterpret = {0};

      reInterpret.fN = number;

      const u64 numberBytes = reInterpret.uN;
            
      const u32 biasExponent = (numberBytes >> 52) & 0x00000000000007FF;
      i32 numberExponent = (i32) biasExponent - 1023;

      const i32 isOdd = numberExponent & 1;
            
      numberExponent -= isOdd;
      
      const u64 normalizedBytes = ((u64) (1023 + isOdd) << 52) | (numberBytes & 0x000FFFFFFFFFFFFFLLU);

      reInterpret.uN = normalizedBytes;

      f64 normalizedNumber = reInterpret.fN;

      f64 try = 0.5 * (1.0 + normalizedNumber);
      
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

f32 SSqrt2F(f32 number, i32 iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      U32F32ToF32U32 reInterpret = {};

      reInterpret.fN = number;

      const u32 numberBytes = reInterpret.uN;
            
      const u32 biasExponent = (numberBytes >> 23) & 0x00000000000000FF;
      i32 numberExponent = (i32) biasExponent - 127;

      const u32 isOdd = numberExponent & 1;
    
      numberExponent -= isOdd;

      const u32 normalizedBytes = ((u32) (127 + isOdd) << 23) | (numberBytes & 0x00000000007FFFFF);

      reInterpret.uN = normalizedBytes;

      f32 normalizedNumber = reInterpret.fN;

      f32 try = 0.5 * (1.0 + normalizedNumber);
      
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

f64 HSqrt(f64 number)
{
      __m128d reg = _mm_set_sd(number);
      __m128d value = _mm_sqrt_sd(reg, reg);

      return _mm_cvtsd_f64(value);
}

f32 HSqrtF(f32 number)
{
      __m128 reg = _mm_set_ss(number);
      __m128 value = _mm_sqrt_ss(reg);

      return _mm_cvtss_f32(value);
}

f32 HFInvSqrtF(f32 number)
{
      __m128 reg = _mm_set_ss(number);
      __m128 value = _mm_rsqrt_ss(reg);

      return _mm_cvtss_f32(value);
}

f32 HFSqrtF(f32 number)
{
      __m128 reg = _mm_set_ss(number);
      __m128 value = _mm_rsqrt_ss(reg);

      return _mm_cvtss_f32(value) * number;
}