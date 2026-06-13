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
      U64F64ToF64U64 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FF0000000000000 + (reInterpret.uN >> 1);

      f64 try = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = 0.5 * (try + (number / try));
      }

      return try;
}

f32 SSqrtF(f32 number, i32 iterations)
{ 
      U32F32ToF32U32 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FBC0000 + (reInterpret.uN >> 1);

      f32 try = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = 0.5f * (try + (number / try));
      }

      return try;
}

f64 SSqrt2(f64 number, i32 iterations)
{ 
      U64F64ToF64U64 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FF0000000000000 + (reInterpret.uN >> 1);

      f64 try = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = try * ((try * try) + (3.0 * number)) / ((3.0 * (try * try)) + number);
      }

      return try;
}

f32 SSqrt2F(f32 number, i32 iterations)
{ 
      U32F32ToF32U32 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FBC0000 + (reInterpret.uN >> 1);

      f32 try = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            try = try * ((try * try) + (3.0f * number)) / ((3.0f * (try * try)) + number);
      }

      return try;
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

f64 SExp(f64 number, f64 exp)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;

      f64 fLogExp = (f64) ((reInterpret.uN >> 52) & 0x7FF) - 1023.0;
      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFFLLU) | 0x3FF0000000000000LLU;

      fLogExp += -2.475404 + reInterpret.fN * (3.321382 + reInterpret.fN * (-0.999602 + 0.153624 * reInterpret.fN));
      fLogExp *= exp;

      const f64 fixedExp = HFloor(fLogExp);

      const i32 intExp = (i32) fixedExp;
      const f64 decExp = fLogExp - fixedExp;

      const f64 decPart = 1.0 + decExp * (0.695713 + decExp * (0.222458 + decExp * 0.081829));
      const u64 finalLogExp = (u64) (intExp + 1023);

      reInterpret.fN = decPart;
      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFFLLU) | (finalLogExp << 52);

      return reInterpret.fN;
}

f32 SExpF(f32 number, f32 exp)
{
      U32F32ToF32U32 reInterpret;
      reInterpret.fN = number;

      f32 fLogExp = (f32) ((reInterpret.uN >> 23) & 0xFF) - 127.0f;
      reInterpret.uN = (reInterpret.uN & 0x007FFFFF) | 0x3F800000;
      
      fLogExp += -2.475404f + reInterpret.fN * (3.321382f + reInterpret.fN * (-0.999602f + 0.153624f * reInterpret.fN));
      fLogExp *= exp;

      const f32 fixedExp = HFloorF(fLogExp);

      const i32 intExp = (i32) fixedExp;
      const f32 decExp = fLogExp - fixedExp;

      const f32 decPart = 1.0f + decExp * (0.695713f + decExp * (0.222458f + decExp * 0.081829f));
      const u32 finalLogExp = (u32) (intExp + 127);

      reInterpret.fN = decPart;
      reInterpret.uN = (reInterpret.uN & 0x007FFFFF) | (finalLogExp << 23);

      return reInterpret.fN;
}

f64 HRound(f64 number)
{
      __m128d reg = _mm_set_sd(number);
      __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

f32 HRoundF(f32 number)
{
      __m128 reg = _mm_set_ss(number);
      __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}

f64 HFloor(f64 number)
{
      __m128d reg = _mm_set_sd(number);
      __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

f32 HFloorF(f32 number)
{
      __m128 reg = _mm_set_ss(number);
      __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}

f64 HCeil(f64 number)
{
      __m128d reg = _mm_set_sd(number);
      __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

f32 HCeilF(f32 number)
{
      __m128 reg = _mm_set_ss(number);
      __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}