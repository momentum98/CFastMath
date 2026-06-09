#include "fastmath.h"

#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>

double SoftwareSqrtFrstLvl(double number, int iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      const u64 numberBytes = *(u64*) &number;
            
      const u32 biasExponent = (numberBytes >> 52) & 0x00000000000007FF;
      i32 numberExponent = (i32) biasExponent - 1023;

      const i32 isOdd = numberExponent & 1;
            
      numberExponent -= isOdd;
            
      const u64 normalizedBytes = ((u64) (1023 + isOdd) << 52) | (numberBytes & 0x000FFFFFFFFFFFFFLLU);
      double normalizedNumber = *(double*) &normalizedBytes;

      double try = 0.5 * (1.0 + normalizedNumber);
      
      for (i32 i = 0; i < iterations; i++)
      {
            try = 0.5 * (try + (normalizedNumber / try));
      }

      const i32 finalExponent = (numberExponent >> 1) + 1023;

      u64 tryBytes = *(u64*) &try;
      u64 finalBytes = ((u64) finalExponent << 52) | (tryBytes & 0x000FFFFFFFFFFFFFLLU);

      return *(double*) &finalBytes;
}

float SoftwareSqrtFrstLvlF(float number, int iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      const u32 numberBytes = *(u32*) &number;
            
      const u32 biasExponent = (numberBytes >> 23) & 0x00000000000000FF;
      i32 numberExponent = (i32) biasExponent - 127;

      const u32 isOdd = numberExponent & 1;
    
      numberExponent -= isOdd;

      const u32 normalizedBytes = ((u32) (127 + isOdd) << 23) | (numberBytes & 0x00000000007FFFFF);
      float normalizedNumber = *(float*) &normalizedBytes;

      float try = 0.5 * (1.0 + normalizedNumber);
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = 0.5 * (try + (normalizedNumber / try));
      }

      const i32 finalExponent = (numberExponent >> 1) + 127;

      u32 tryBytes = *(u32*) &try;
      u32 finalBytes = ((u32) finalExponent << 23) | (tryBytes & 0x00000000007FFFFF);

      return *(float*) &finalBytes;
}

double SoftwareSqrtScndLvl(double number, int iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      const u64 numberBytes = *(u64*) &number;
            
      const u32 biasExponent = (numberBytes >> 52) & 0x00000000000007FF;
      i32 numberExponent = (i32) biasExponent - 1023;

      const i32 isOdd = numberExponent & 1;
    
      numberExponent -= isOdd;

      const u64 normalizedBytes = ((u64) (1023 + isOdd) << 52) | (numberBytes & 0x000FFFFFFFFFFFFFLLU);
      double normalizedNumber = *(double*) &normalizedBytes;

      double try = 0.5 * (1.0 + normalizedNumber);
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = try * ((try * try) + (3 * normalizedNumber)) / ((3 * (try * try)) + normalizedNumber);
      }

      const i32 finalExponent = (numberExponent >> 1) + 1023;

      u64 tryBytes = *(u64*) &try;
      u64 finalBytes = ((u64) finalExponent << 52) | (tryBytes & 0x000FFFFFFFFFFFFFLLU);

      return *(double*) &finalBytes;
}

float SoftwareSqrtScndLvlF(float number, int iterations)
{ 
      if (number <= 0)
      {
            return number;
      }

      const u32 numberBytes = *(u32*) &number;
            
      const u32 biasExponent = (numberBytes >> 23) & 0x00000000000000FF;
      i32 numberExponent = (i32) biasExponent - 127;

      const u32 isOdd = numberExponent & 1;
    
      numberExponent -= isOdd;

      const u32 normalizedBytes = ((u32) (127 + isOdd) << 23) | (numberBytes & 0x00000000007FFFFF);
      float normalizedNumber = *(float*) &normalizedBytes;

      float try = 0.5 * (1.0 + normalizedNumber);
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = try * ((try * try) + (3 * normalizedNumber)) / ((3 * (try * try)) + normalizedNumber);
      }

      const i32 finalExponent = (numberExponent >> 1) + 127;

      u32 tryBytes = *(u32*) &try;
      u32 finalBytes = ((u32) finalExponent << 23) | (tryBytes & 0x00000000007FFFFF);

      return *(float*) &finalBytes;
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

      return (_mm_cvtss_f32(value)) * number;
}