#include "fastmath.h"

#include <immintrin.h>

const u32 OCTANT_FLAGS = 0x66;

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

f64 HRound(f64 number)
{
      const __m128d reg = _mm_set_sd(number);
      const __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

f32 HRoundF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}

f64 HFloor(f64 number)
{
      const __m128d reg = _mm_set_sd(number);
      const __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

f32 HFloorF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}

f64 HCeil(f64 number)
{
      const __m128d reg = _mm_set_sd(number);
      const __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

f32 HCeilF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}

f64 SAbs(f64 number)
{
      U64F64ToF64U64 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = (reInterpret.uN & 0x7FFFFFFFFFFFFFFF);

      return reInterpret.fN;
}

f32 SAbsF(f32 number)
{
      U32F32ToF32U32 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = (reInterpret.uN & 0x7FFFFFFF);

      return reInterpret.fN;
}

f64 SCopySign(f64 number, f64 signNumber)
{
      U64F64ToF64U64 reInterpretNum;
      U64F64ToF64U64 reInterpretSign;

      reInterpretSign.fN = signNumber;
      reInterpretNum.fN = number;
      
      reInterpretNum.uN = (reInterpretNum.uN & 0x7FFFFFFFFFFFFFFFLLU) | (reInterpretSign.uN & 0x8000000000000000LLU);

      return reInterpretNum.fN;
}

f32 SCopySignF(f32 number, f32 signNumber)
{
      U32F32ToF32U32 reInterpretNum;
      U32F32ToF32U32 reInterpretSign;

      reInterpretSign.fN = signNumber;
      reInterpretNum.fN = number;

      reInterpretNum.uN = (reInterpretNum.uN & 0x7FFFFFFF) | (reInterpretSign.uN & 0x80000000);

      return reInterpretNum.fN;
}

f64 SInv(f64 number, u32 iterations)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;

      reInterpret.uN = 0x7FDE648000000000LLU - reInterpret.uN; 
      f64 fValue = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            fValue = fValue * (2.0 - number * fValue);
      }

      return fValue;
}

f32 SInvF(f32 number, u32 iterations)
{
      U32F32ToF32U32 reInterpret;
      reInterpret.fN = number;

      reInterpret.uN = 0x7EF311C2 - reInterpret.uN; 
      f32 fValue = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            fValue = fValue * (2.0f - number * fValue);
      }

      return fValue;
}

f64 HInv(f64 number, u32 iterations)
{
      const __m128 reg = _mm_set_ss((f32) number);
      const __m128 value = _mm_rcp_ss(reg);

      f64 fValue = (f64) _mm_cvtss_f32(value);

      for (u32 i = 0; i < iterations; ++i)
      {
            fValue = fValue * (2.0 - number * fValue);
      }

      return fValue;
}

f32 HInvF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_rcp_ss(reg);

      return _mm_cvtss_f32(value);
}

f64 SSqrt(f64 number, i32 iterations)
{ 
      U64F64ToF64U64 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FF0000000000000 + (reInterpret.uN >> 1);

      f64 guess = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            guess = 0.5 * (guess + (number * SInv(guess, 3)));
      }

      return guess;
}

f32 SSqrtF(f32 number, i32 iterations)
{ 
      U32F32ToF32U32 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FBC0000 + (reInterpret.uN >> 1);

      f32 guess = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            guess = 0.5f * (guess + (number * SInvF(guess, 3)));
      }

      return guess;
}

f64 SSqrt2(f64 number, i32 iterations)
{ 
      U64F64ToF64U64 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FF0000000000000 + (reInterpret.uN >> 1);

      f64 guess = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            guess = guess * ((guess * guess) + (3.0 * number)) * SInv((3.0 * (guess * guess)) + number, 3);
      }

      return guess;
}

f32 SSqrt2F(f32 number, i32 iterations)
{ 
      U32F32ToF32U32 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FBC0000 + (reInterpret.uN >> 1);

      f32 guess = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            guess = guess * ((guess * guess) + (3.0f * number)) * SInvF((3.0f * (guess * guess)) + number, 3);
      }

      return guess;
}

f64 HSqrt(f64 number)
{
      const __m128d reg = _mm_set_sd(number);
      const __m128d value = _mm_sqrt_sd(reg, reg);

      return _mm_cvtsd_f64(value);
}

f32 HSqrtF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_sqrt_ss(reg);

      return _mm_cvtss_f32(value);
}

f32 HFInvSqrtF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_rsqrt_ss(reg);

      return _mm_cvtss_f32(value);
}

f32 HFSqrtF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_rsqrt_ss(reg);

      return _mm_cvtss_f32(value) * number;
}

f64 SExp(f64 number, f64 exp)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;

      f64 fLogExp = (f64) ((reInterpret.uN >> 52) & 0x7FF) - 1023.0;
      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFFLLU) | 0x3FF0000000000000LLU;

      fLogExp += -3.0390146 + reInterpret.fN * (4.9883594 + reInterpret.fN * (-3.1583092 + reInterpret.fN * (1.4925010 + reInterpret.fN * (-0.4211116 + 0.0515152 * reInterpret.fN))));
      fLogExp *= exp;

      const f64 fixedExp = HFloor(fLogExp);

      const i32 intExp = (i32) fixedExp;
      const f64 decExp = fLogExp - fixedExp;

      const f64 decPart = 1.0 + decExp * (0.69314718 + decExp * (0.24022647 + decExp * (0.05550513 + decExp * (0.00961402 + 0.00134234 * decExp))));
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
      
      fLogExp += -3.0390146f + reInterpret.fN * (4.9883594f + reInterpret.fN * (-3.1583092f + reInterpret.fN * (1.4925010f + reInterpret.fN * (-0.4211116f + 0.0515152f * reInterpret.fN))));
      fLogExp *= exp;

      const f32 fixedExp = HFloorF(fLogExp);

      const i32 intExp = (i32) fixedExp;
      const f32 decExp = fLogExp - fixedExp;

      const f32 decPart = 1.0f + decExp * (0.69314718f + decExp * (0.24022647f + decExp * (0.05550513f + decExp * (0.00961402f + 0.00134234f * decExp))));
      const u32 finalLogExp = (u32) (intExp + 127);

      reInterpret.fN = decPart;
      reInterpret.uN = (reInterpret.uN & 0x007FFFFF) | (finalLogExp << 23);

      return reInterpret.fN;
}

f64 SCos(f64 turn)
{
      const f64 interval = (turn + 0.25) - HFloorF(turn + 0.25);
      const f64 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f64 fIsOdd = (f64) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 2) & 1);

      const f64 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1;
      const f64 fUseCos = (f64) useCos;
      
      const f64 cInterval4 = cIntervalSqr * cIntervalSqr;

      const f64 res = (useCos) ? (1.0 - cIntervalSqr * 0.308425137 + cInterval4 * (0.015854344 - 0.000325991 * cIntervalSqr)) : (cInterval * (0.785398163 - cIntervalSqr * 0.080745512 + cInterval4 * (0.002490095 - 0.000036584 * cIntervalSqr)));
      
      return res * fSign;
}

f32 SCosF(f32 turn)
{
      const f32 interval = (turn + 0.25f) - HFloorF(turn + 0.25f);
      const f32 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f32 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f32 fIsOdd = (f32) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0f - 2.0f * fIsOdd));

      const f32 fSign = 1.0f - 2.0f * (f32) ((octant >> 2) & 1);

      const f32 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1; 
      const f32 fUseCos = (f32) useCos;
      
      const f32 cInterval4 = cIntervalSqr * cIntervalSqr;

      const f32 res = (useCos) ? (1.0f - cIntervalSqr * 0.308425137f + cInterval4 * (0.015854344f - 0.000325991f * cIntervalSqr)) : (cInterval * (0.785398163f - cIntervalSqr * 0.080745512f + cInterval4 * (0.002490095f - 0.000036584f * cIntervalSqr)));
      
      return res * fSign;
}

f64 SSin(f64 turn)
{
      const f64 interval = turn - HFloorF(turn);
      const f64 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f64 fIsOdd = (f64) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 2) & 1);

      const f64 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1;
      const f64 fUseCos = (f64) useCos;
      
      const f64 cInterval4 = cIntervalSqr * cIntervalSqr;

      const f64 res = (useCos) ? (1.0 - cIntervalSqr * 0.308425137 + cInterval4 * (0.015854344 - 0.000325991 * cIntervalSqr)) : (cInterval * (0.785398163 - cIntervalSqr * 0.080745512 + cInterval4 * (0.002490095 - 0.000036584 * cIntervalSqr)));
      
      return res * fSign;
}

f32 SSinF(f32 turn)
{
      const f32 interval = turn - HFloorF(turn);
      const f32 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f32 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f32 fIsOdd = (f32) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0f - 2.0f * fIsOdd));

      const f32 fSign = 1.0f - 2.0f * (f32) ((octant >> 2) & 1);

      const f32 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1; 
      const f32 fUseCos = (f32) useCos;

      const f32 cInterval4 = cIntervalSqr * cIntervalSqr;

      const f32 res = (useCos) ? (1.0f - cIntervalSqr * 0.308425137f + cInterval4 * (0.015854344f - 0.000325991f * cIntervalSqr)) : (cInterval * (0.785398163f - cIntervalSqr * 0.080745512f + cInterval4 * (0.002490095f - 0.000036584f * cIntervalSqr)));
      
      return res * fSign;
}

f64 SFCos(f64 turn)
{
      const f64 interval = (turn + 0.25) - HFloorF(turn + 0.25);
      const f64 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f64 fIsOdd = (f64) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 2) & 1);

      f64 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1;
      const f64 fUseCos = (f64) useCos;

      const f64 res = (useCos) ? (1.0 - 0.30842 * cIntervalSqr) : (cInterval * (0.78539 - 0.08074 * cIntervalSqr));
      
      return res * fSign;
}

f32 SFCosF(f32 turn)
{
      const f32 interval = (turn + 0.25f) - HFloorF(turn + 0.25f);
      const f32 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f32 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f32 fIsOdd = (f32) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0f - 2.0f * fIsOdd));

      const f32 fSign = 1.0f - 2.0f * (f32) ((octant >> 2) & 1);

      f32 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1; 
      const f32 fUseCos = (f32) useCos;

      const f32 res = (useCos) ? (1.0f - 0.30842f * cIntervalSqr) : (cInterval * (0.78539f - 0.08074f * cIntervalSqr));
      
      return res * fSign;
}

f64 SFSin(f64 turn)
{
      const f64 interval = turn - HFloorF(turn);
      const f64 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f64 fIsOdd = (f64) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 2) & 1);

      f64 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1;
      const f64 fUseCos = (f64) useCos;

      const f64 res = (useCos) ? (1.0 - 0.30842 * cIntervalSqr) : (cInterval * (0.78539 - 0.08074 * cIntervalSqr));
      
      return res * fSign;
}

f32 SFSinF(f32 turn)
{
      const f32 interval = turn - HFloorF(turn);
      const f32 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f32 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f32 fIsOdd = (f32) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0f - 2.0f * fIsOdd));

      const f32 fSign = 1.0f - 2.0f * (f32) ((octant >> 2) & 1);

      f32 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1; 
      const f32 fUseCos = (f32) useCos;
      
      const f32 res = (useCos) ? (1.0f - 0.30842f * cIntervalSqr) : (cInterval * (0.78539f - 0.08074f * cIntervalSqr));

      return res * fSign;
}

f64 SFTan(f64 turn)
{
      const f64 interval = turn - HFloorF(turn);
      const f64 qInterval = interval * 8.0;
      const u32 octant = (u32) HFloorF(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;
      const f64 fIsOdd = (f64) isOdd;

      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      f64 cIntervalSqr = cInterval * cInterval;

      f64 result = cInterval * (0.785398163 + cIntervalSqr * (0.16137528 + 0.05322692 * cIntervalSqr));

      result = isOdd ? HInv(result, 1) : result;

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 1) & 1);

      return result * fSign;
}

f32 SFTanF(f32 turn)
{
      const f32 interval = turn - HFloorF(turn);
      const f32 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f32 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;
      const f32 fIsOdd = (f32) isOdd;

      cInterval = fIsOdd + (cInterval * (1.0f - 2.0f * fIsOdd));

      f32 cIntervalSqr = cInterval * cInterval;

      f32 result = cInterval * (0.785398163f + cIntervalSqr * (0.16137528f + 0.05322692f * cIntervalSqr));

      result = isOdd ? HInvF(result) : result;

      const f32 fSign = 1.0f - 2.0f * (f32) ((octant >> 1) & 1);

      return result * fSign;
}

f64 STan(f64 turn)
{
      const f64 interval = turn - HFloorF(turn);
      const f64 qInterval = interval * 8.0;
      const u32 octant = (u32) HFloorF(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;
      const f64 fIsOdd = (f64) isOdd;

      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      f64 cIntervalSqr = cInterval * cInterval;

      f64 result = cInterval * (0.785398163 + cIntervalSqr * (0.16137528 + 0.05322692 * cIntervalSqr));

      result = isOdd ? (1.0 / result) : result;

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 1) & 1);

      return result * fSign;
}

f32 STanF(f32 turn)
{
      const f32 interval = turn - HFloorF(turn);
      const f32 qInterval = interval * 8.0f;
      const u32 octant = (u32) HFloorF(qInterval);

      f32 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;
      const f32 fIsOdd = (f32) isOdd;

      cInterval = fIsOdd + (cInterval * (1.0f - 2.0f * fIsOdd));

      f32 cIntervalSqr = cInterval * cInterval;

      f32 result = cInterval * (0.785398163f + cIntervalSqr * (0.16137528f + 0.05322692f * cIntervalSqr));

      result = isOdd ? (1.0f / result) : result;

      const f32 fSign = 1.0f - 2.0f * (f32) ((octant >> 1) & 1);

      return result * fSign;
}

f64 SAtan2(f64 x, f64 y)
{
      const f64 absX = SAbs(x);
      const f64 absY = SAbs(y);
      
      const u32 cond = absY > absX;

      const f64 a = cond ? absX : absY;
      const f64 b = cond ? absY : absX;

      const f64 slope = a * SInv(b, 2);
      
      f64 slopeSqr = slope * slope;
      f64 angle = (slopeSqr * -0.19194795 + 0.97239411) * slope;
      
      angle = cond ? 1.57079632 - angle : angle;
      angle = (x < 0.0) ? 3.14159265 - angle : angle;
      
      return SCopySign(angle, y); 
}

f32 SAtan2F(f32 x, f32 y)
{
      const f32 absX = SAbsF(x);
      const f32 absY = SAbsF(y);
      
      const u32 cond = absY > absX;

      const f32 a = cond ? absX : absY;
      const f32 b = cond ? absY : absX;

      const f32 slope = a * SInvF(b, 2);
      
      f32 slopeSqr = slope * slope;
      f32 angle = (slopeSqr * -0.19194795f + 0.97239411f) * slope;
      
      angle = cond ? 1.57079632f - angle : angle;
      angle = (x < 0.0f) ? 3.14159265f - angle : angle;
      
      return SCopySign(angle, y); 
}