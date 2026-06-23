#include "fastmath.h"

#include <immintrin.h>

const u32 OCTANT_FLAGS = 0x66;

const f64 PI = 3.14159265358979323846264338327950288;
const f32 PIF = 3.14159265358979323846f;

const f64 INV_PI = 0.31830988618379067154;
const f32 INV_PIF = 0.31830988618f;

const f64 INV_FULL_ANGLE = 1.0 / 360.0;
const f32 INV_FULL_ANGLEF = 1.0f / 360.0f;

//
// TODO: Change All Small Methods Type to Inline.
//

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
      f32 fN;
      i32 uN;
} I32F32ToF32I32;

typedef union
{
      f64 fN;
      i64 uN;
} I64F64ToF64I64;

f64 WrapTo180(f64 angle)
{
      return angle - 360.0 * (HFloor((angle + 180.0) * INV_FULL_ANGLE));
}

f32 WrapTo180F(f32 angle)
{
      return angle - (360.0f * HFloorF((angle + 180.0f) * INV_FULL_ANGLEF));
}

f64 WrapTo360(f64 angle)
{
      return angle - (360.0 * HFloor(angle * INV_FULL_ANGLE));
}

f32 WrapTo360F(f32 angle)
{
      return angle - (360.0f * HFloorF(angle * INV_FULL_ANGLEF));
}

f64 RadToDeg(f64 rad)
{
      return (rad * INV_PI) * 180.0;
}

f32 RadToDegF(f32 rad)
{
      return (rad * INV_PIF) * 180.0f;
}

f64 RadToTurn(f64 rad)
{
      return (rad * INV_PI * 0.5);
}

f32 RadToTurnF(f32 rad)
{
      return (rad * INV_PIF * 0.5f);
}

f64 DegToRad(f64 deg)
{
      return (deg * INV_FULL_ANGLE) * (2.0 * PI);
}

f32 DegToRadF(f32 deg)
{
      return (deg * INV_FULL_ANGLEF) * (2.0f * PIF);
}

f64 DegToTurn(f64 deg)
{
      return (deg * INV_FULL_ANGLE);
}

f32 DegToTurnF(f32 deg)
{
      return (deg * INV_FULL_ANGLEF);
}

f64 TurnToRad(f64 turn)
{
      return turn * (2.0 * PI);
}

f32 TurnToRadF(f32 turn)
{
      return turn * (2.0f * PIF);
}

f64 TurnToDeg(f64 turn)
{
      return turn * 360.0;
}

f32 TurnToDegF(f32 turn)
{
      return turn * 360.0f;
}

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

u32 IsNegative(f64 number)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;
      
      return (reInterpret.uN >> 63) & 1;
}

u32 IsNegativeF(f32 number)
{
      U32F32ToF32U32 reInterpret;
      reInterpret.fN = number;
      
      return (reInterpret.uN >> 31) & 1;
}

i32 SSign(f64 number)
{
      return (number > 0.0) - (number < 0.0);
}

i32 SSignF(f32 number)
{
      return (number > 0.0f) - (number < 0.0f);
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
      
      reInterpretNum.uN = (reInterpretNum.uN & 0x7FFFFFFFFFFFFFFF) | (reInterpretSign.uN & 0x8000000000000000);

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

f64 SFLerp(f64 a, f64 b, f64 p)
{
      return a + (b - a) * p;
}

f32 SFLerpF(f32 a, f32 b, f32 p)
{
      return a + (b - a) * p;
}

f64 SLerp(f64 a, f64 b, f64 p)
{
      return (1 - p) * a + b * p;
}

f32 SLerpF(f32 a, f32 b, f32 p)
{
      return (1 - p) * a + b * p;
}

f64 SEaseIn(f64 a, f64 b, f64 p)
{
      const f64 finalProgress = (b - a) * (p * p);
      
      return a + finalProgress;
}

f32 SEaseInF(f32 a, f32 b, f32 p)
{
      const f32 finalProgress = (b - a) * (p * p);
      
      return a + finalProgress;
}

f64 SEaseOut(f64 a, f64 b, f64 p)
{
      const f64 finalProgress = (b - a) * (p * (2.0 - p));
      
      return a + finalProgress;
}

f32 SEaseOutF(f32 a, f32 b, f32 p)
{
      const f32 finalProgress = (b - a) * (p * (2.0f - p));
      
      return a + finalProgress;
}

f64 SEaseInOut(f64 a, f64 b, f64 p)
{
      const f64 finalProgress = (b - a) * (p * p * (3.0 - 2.0 * p));
      
      return a + finalProgress;
}

f32 SEaseInOutF(f32 a, f32 b, f32 p)
{
      const f32 finalProgress = (b - a) * (p * p * (3.0f - 2.0f * p));
      
      return a + finalProgress;
}

f64 SInv(f64 number, u32 iterations)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;

      reInterpret.uN = 0x7FDE648000000000 - reInterpret.uN; 
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

f32 HInvF(f32 number, u32 iterations)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_rcp_ss(reg);

      f32 fValue = _mm_cvtss_f32(value);

      for (u32 i = 0; i < iterations; ++i)
      {
            fValue = fValue * (2.0f - number * fValue);
      }

      return fValue;
}

f64 SSqrt(f64 number, u32 iterations)
{ 
      U64F64ToF64U64 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FF0000000000000 + (reInterpret.uN >> 1);

      f64 guess = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            guess = 0.5 * (guess + (number / guess));
      }

      return guess;
}

f32 SSqrtF(f32 number, u32 iterations)
{ 
      U32F32ToF32U32 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FBC0000 + (reInterpret.uN >> 1);

      f32 guess = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            guess = 0.5f * (guess + (number / guess));
      }

      return guess;
}

f64 SSqrt2(f64 number, u32 iterations)
{ 
      U64F64ToF64U64 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FF0000000000000 + (reInterpret.uN >> 1);

      f64 guess = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            guess = guess * ((guess * guess) + (3.0 * number)) / (3.0 * (guess * guess));
      }

      return guess;
}

f32 SSqrt2F(f32 number, u32 iterations)
{ 
      U32F32ToF32U32 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = 0x1FBC0000 + (reInterpret.uN >> 1);

      f32 guess = reInterpret.fN;
      
      for (u32 i = 0; i < iterations; ++i)
      {
            guess = guess * ((guess * guess) + (3.0f * number)) / (3.0 * (guess * guess));
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
      
      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFF) | 0x3FF0000000000000;
      
      const f64 nMantissa = reInterpret.fN - 1.0;

      fLogExp += nMantissa * (1.4426950408889634074 + nMantissa * (
            -0.7213475153215513511 + nMantissa * (
             0.4808983455799981395 + nMantissa * (
            -0.3606740620216656730 + nMantissa * (
             0.2885437812959828591 + nMantissa * (
            -0.2401140026262445199 + nMantissa * (
             0.2061214040182851410 + nMantissa * (
            -0.1774395687796900400 + nMantissa * (
             0.1545620942548261300 + nMantissa * (
            -0.1287968412850722200 + nMantissa * 0.1068478440700465500
      ))))))))));

      fLogExp *= exp;

      const f64 fixedExp = HFloor(fLogExp);
      const i32 intExp = (i32) fixedExp;
      const f64 decExp = fLogExp - fixedExp;

      const f64 decPart = 1.0 + decExp * (
            0.6931471805599453094 + decExp * (
            0.2402265069591007183 + decExp * (
            0.0555041086648215799 + decExp * (
            0.0096181291076284776 + decExp * (
            0.0013333558146428443 + decExp * (
            0.0001540353039338160 + decExp * (
            0.0000152527338048598 + decExp * (
            0.0000013215486790144 + decExp * 0.0000001017808600924
      ))))))));

      const u64 finalLogExp = (u64) (intExp + 1023);

      reInterpret.fN = decPart;
      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFF) | (finalLogExp << 52);

      return reInterpret.fN;
}

f32 SExpF(f32 number, f32 exp)
{
      U32F32ToF32U32 reInterpret;
      reInterpret.fN = number;

      f32 fLogExp = (f32) ((reInterpret.uN >> 23) & 0xFF) - 127.0f;
      reInterpret.uN = (reInterpret.uN & 0x007FFFFF) | 0x3F800000;
      
      const f32 nMantissa = reInterpret.fN - 1.0f;

      fLogExp += nMantissa * (1.4426950f + nMantissa * (
            -0.7211625f + nMantissa * (
             0.4731114f + nMantissa * (
            -0.3150020f + nMantissa * (
             0.1987642f + nMantissa * (
            -0.1171452f + nMantissa * 0.0524026f
      ))))));

      fLogExp *= exp;

      const f32 fixedExp = HFloorF(fLogExp);
      const i32 intExp = (i32) fixedExp;
      const f32 decExp = fLogExp - fixedExp;

      const f32 decPart = 1.0f + decExp * (
            0.6931472f + decExp * (
            0.2402265f + decExp * (
            0.0555041f + decExp * (
            0.0096181f + decExp * (
            0.0013334f + decExp * (
            0.0001540f + decExp * 0.0000153f
      ))))));

      const u32 finalLogExp = (u32) (intExp + 127);

      reInterpret.fN = decPart;
      reInterpret.uN = (reInterpret.uN & 0x007FFFFF) | (finalLogExp << 23);

      return reInterpret.fN;
}

f64 SMax(f64 a, f64 b)
{
      return ((a + b) + SAbs(a - b)) * 0.5;
}

f32 SMaxF(f32 a, f32 b)
{
      return ((a + b) + SAbsF(a - b)) * 0.5f;
}

f64 HMax(f64 a, f64 b)
{
      const __m128d aReg = _mm_set_sd(a);
      const __m128d bReg = _mm_set_sd(b);

      const __m128d value = _mm_max_sd(aReg, bReg);

      return _mm_cvtsd_f64(value);
}

f32 HMaxF(f32 a, f32 b)
{
      const __m128 aReg = _mm_set_ss(a);
      const __m128 bReg = _mm_set_ss(b);

      const __m128 value = _mm_max_ss(aReg, bReg);

      return _mm_cvtss_f32(value);
}

f64 SMin(f64 a, f64 b)
{
      return ((a + b) - SAbs(a - b)) * 0.5;
}

f32 SMinF(f32 a, f32 b)
{
      return ((a + b) - SAbsF(a - b)) * 0.5f;
}

f64 HMin(f64 a, f64 b)
{
      const __m128d aReg = _mm_set_sd(a);
      const __m128d bReg = _mm_set_sd(b);

      const __m128d value = _mm_min_sd(aReg, bReg);

      return _mm_cvtsd_f64(value);
}

f32 HMinF(f32 a, f32 b)
{
      const __m128 aReg = _mm_set_ss(a);
      const __m128 bReg = _mm_set_ss(b);

      const __m128 value = _mm_min_ss(aReg, bReg);

      return _mm_cvtss_f32(value);
}

f64 SClamp(f64 value, f64 min, f64 max)
{
      const f64 minPatched = SMax(value, min);
      const f64 maxPatched = SMin(minPatched, max);

      return maxPatched;
}

f32 SClampF(f32 value, f32 min, f32 max)
{
      const f32 minPatched = SMaxF(value, min);
      const f32 maxPatched = SMinF(minPatched, max);

      return maxPatched;
}

f64 HClamp(f64 value, f64 min, f64 max)
{
      const f64 minPatched = HMax(value, min);
      const f64 maxPatched = HMin(minPatched, max);

      return maxPatched;
}

f32 HClampF(f32 value, f32 min, f32 max)
{
      const f32 minPatched = HMaxF(value, min);
      const f32 maxPatched = HMinF(minPatched, max);

      return maxPatched;
}

f64 SFExp(f64 number, f64 exp)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;

      f64 fLogExp = (f64) ((reInterpret.uN >> 52) & 0x7FF) - 1023.0;
      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFF) | 0x3FF0000000000000;

      fLogExp += -3.0390146 + reInterpret.fN * (4.9883594 + reInterpret.fN * (-3.1583092 + reInterpret.fN * (1.4925010 + reInterpret.fN * (-0.4211116 + 0.0515152 * reInterpret.fN))));
      fLogExp *= exp;

      const f64 fixedExp = HFloor(fLogExp);

      const i32 intExp = (i32) fixedExp;
      const f64 decExp = fLogExp - fixedExp;

      const f64 decPart = 1.0 + decExp * (0.69314718 + decExp * (0.24022647 + decExp * (0.05550513 + decExp * (0.00961402 + 0.00134234 * decExp))));
      const u64 finalLogExp = (u64) (intExp + 1023);

      reInterpret.fN = decPart;
      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFF) | (finalLogExp << 52);

      return reInterpret.fN;
}

f32 SFExpF(f32 number, f32 exp)
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
      const f64 interval = (turn + 0.25) - HFloor(turn + 0.25);
      const f64 qInterval = interval * 8.0;
      const u32 octant = (u32) HFloor(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f64 fIsOdd = (f64) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 2) & 1);

      const f64 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1;
      const f64 fUseCos = (f64) useCos;
      
      const f64 cInterval4 = cIntervalSqr * cIntervalSqr;

      const f64 cosPart = 1.0 - cIntervalSqr * 0.308425137 + cInterval4 * (0.015854344 - 0.000325991 * cIntervalSqr);
      const f64 sinPart = cInterval * (0.785398163 - cIntervalSqr * 0.080745512 + cInterval4 * (0.002490095 - 0.000036584 * cIntervalSqr));

      const f64 res = useCos * cosPart + (f64) (1 - useCos) * sinPart;

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
      
      const f32 cInterval4 = cIntervalSqr * cIntervalSqr;

      const f32 cosPart = 1.0f - cIntervalSqr * 0.308425137f + cInterval4 * (0.015854344f - 0.000325991f * cIntervalSqr);
      const f32 sinPart = cInterval * (0.785398163f - cIntervalSqr * 0.080745512f + cInterval4 * (0.002490095f - 0.000036584f * cIntervalSqr));

      const f32 res = useCos * cosPart + (f32) (1 - useCos) * sinPart;
      
      return res * fSign;
}

f64 SSin(f64 turn)
{
      const f64 interval = turn - HFloor(turn);
      const f64 qInterval = interval * 8.0;
      const u32 octant = (u32) HFloor(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f64 fIsOdd = (f64) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 2) & 1);

      const f64 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1;
      
      const f64 cInterval4 = cIntervalSqr * cIntervalSqr;

      const f64 cosPart = 1.0 - cIntervalSqr * 0.308425137 + cInterval4 * (0.015854344 - 0.000325991 * cIntervalSqr);
      const f64 sinPart = cInterval * (0.785398163 - cIntervalSqr * 0.080745512 + cInterval4 * (0.002490095 - 0.000036584 * cIntervalSqr));

      const f64 res = useCos * cosPart + (f64) (1 - useCos) * sinPart;

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
      
      const f32 cInterval4 = cIntervalSqr * cIntervalSqr;

      const f32 sinPart = cInterval * (0.785398163f - cIntervalSqr * 0.080745512f + cInterval4 * (0.002490095f - 0.000036584f * cIntervalSqr));
      const f32 cosPart = 1.0f - cIntervalSqr * 0.308425137f + cInterval4 * (0.015854344f - 0.000325991f * cIntervalSqr);

      const f32 res = useCos * cosPart + (f32) (1 - useCos) * sinPart;

      return res * fSign;
}

f64 SFCos(f64 turn)
{
      const f64 interval = (turn + 0.25) - HFloor(turn + 0.25);
      const f64 qInterval = interval * 8.0;
      const u32 octant = (u32) HFloor(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f64 fIsOdd = (f64) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 2) & 1);

      f64 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1;
      
      const f64 cosPart = 1.0 - 0.30842 * cIntervalSqr;
      const f64 sinPart = cInterval * (0.78539 - 0.08074 * cIntervalSqr);

      const f64 res = useCos * cosPart + (f64) (1 - useCos) * sinPart;
      
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
      
      const f32 cosPart = 1.0f - 0.30842f * cIntervalSqr;
      const f32 sinPart = cInterval * (0.78539f - 0.08074f * cIntervalSqr);

      const f32 res = useCos * cosPart + (f32) (1 - useCos) * sinPart;
      
      return res * fSign;
}

f64 SFSin(f64 turn)
{
      const f64 interval = turn - HFloor(turn);
      const f64 qInterval = interval * 8.0;
      const u32 octant = (u32) HFloor(qInterval);

      f64 cInterval = qInterval - octant;

      const u32 isOdd = octant & 1;

      const f64 fIsOdd = (f64) isOdd;
      cInterval = fIsOdd + (cInterval * (1.0 - 2.0 * fIsOdd));

      const f64 fSign = 1.0 - 2.0 * (f64) ((octant >> 2) & 1);

      f64 cIntervalSqr = cInterval * cInterval;

      const u32 useCos = (OCTANT_FLAGS >> octant) & 1;
      
      const f64 cosPart = 1.0 - 0.30842 * cIntervalSqr;
      const f64 sinPart = cInterval * (0.78539 - 0.08074 * cIntervalSqr);

      const f64 res = useCos * cosPart + (f64) (1 - useCos) * sinPart;
      
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
      
      const f32 cosPart = 1.0f - 0.30842f * cIntervalSqr;
      const f32 sinPart = cInterval * (0.78539f - 0.08074f * cIntervalSqr);
      
      const f32 res = useCos * cosPart + (f32) (1 - useCos) * sinPart;

      return res * fSign;
}

f64 SFTan(f64 turn)
{
      const f64 interval = turn - HFloor(turn);
      const f64 qInterval = interval * 8.0;
      const u32 octant = (u32) HFloor(qInterval);

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

      result = isOdd ? HInvF(result, 2) : result;

      const f32 fSign = 1.0f - 2.0f * (f32) ((octant >> 1) & 1);

      return result * fSign;
}

f64 STan(f64 turn)
{
      const f64 interval = turn - HFloor(turn);
      const f64 qInterval = interval * 8.0;
      const u32 octant = (u32) HFloor(qInterval);

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
      
      return SCopySignF(angle, y); 
}

f64 SACos(f64 cos)
{
      const u32 signal = SSign(cos);
      cos = SAbs(cos);

      const f64 cosSqr = cos * cos;
      const f64 cosSqr2 = (cosSqr * cosSqr);
      const f64 cosSqr3 = (cosSqr2 * cosSqr);
      const f64 cosSqr4 = (cosSqr3 * cosSqr);
      const f64 cosSqr5 = (cosSqr4 * cosSqr);

      const f64 p1 = (-0.045423851531649318 * cosSqr4) 
         + (0.74990235515610817 * cosSqr3) 
         - (3.2468404245995405 * cosSqr2) 
         + (5.4637449576887556 * cosSqr) 
         - 3.2556571386590153;

      const f64 p2 = cosSqr5
         - (1.3451543781223945 * cosSqr4) 
         + (8.92015291242372 * cosSqr3) 
         - (26.733221538350166 * cosSqr2) 
         + (37.051740033061225 * cosSqr) 
         - 19.533942831954092;

      return (PI * 0.5) - signal * ((cos + (cos * cosSqr * (p1 * (HInv(p2, 2))))));
}

f32 SACosF(f32 cos)
{
      const u32 signal = SSignF(cos);
      cos = SAbsF(cos);
      
      const f32 cosSqr = cos * cos;
      const f32 cosSqr2 = (cosSqr * cosSqr);
      const f32 cosSqr3 = (cosSqr2 * cosSqr);
      const f32 cosSqr4 = (cosSqr3 * cosSqr);
      const f32 cosSqr5 = (cosSqr4 * cosSqr);

      const f32 p1 = (-0.045423f * cosSqr4) 
         + (0.74990f * cosSqr3) 
         - (3.24684f * cosSqr2) 
         + (5.46374f * cosSqr) 
         - 3.25565f;

      const f32 p2 = cosSqr5
         - (1.34515f * cosSqr4) 
         + (8.92015f * cosSqr3) 
         - (26.73322f * cosSqr2) 
         + (37.05174f * cosSqr)
         - 19.53394f;

      return (PIF * 0.5f) - signal * ((cos + (cos * cosSqr * (p1 * (HInvF(p2, 2))))));
}

f64 SASin(f64 sin)
{
      const u32 signal = SSign(sin);
      sin = SAbs(sin);
      
      const f64 sinSqr = sin * sin;
      const f64 sinSqr2 = (sinSqr * sinSqr);
      const f64 sinSqr3 = (sinSqr2 * sinSqr);
      const f64 sinSqr4 = (sinSqr3 * sinSqr);
      const f64 sinSqr5 = (sinSqr4 * sinSqr);

      const f64 p1 = (-0.045423851531649318 * sinSqr4) 
         + (0.74990235515610817 * sinSqr3) 
         - (3.2468404245995405 * sinSqr2) 
         + (5.4637449576887556 * sinSqr) 
         - 3.2556571386590153;

      const f64 p2 = sinSqr5
         - (1.3451543781223945 * sinSqr4) 
         + (8.92015291242372 * sinSqr3) 
         - (26.733221538350166 * sinSqr2) 
         + (37.051740033061225 * sinSqr) 
         - 19.533942831954092;

      return signal * (sin + (sin * sinSqr * (p1 * (HInv(p2, 2)))));
}

f32 SASinF(f32 sin)
{
      const u32 signal = SSignF(sin);
      sin = SAbsF(sin);
      
      const f32 sinSqr = sin * sin;
      const f32 sinSqr2 = (sinSqr * sinSqr);
      const f32 sinSqr3 = (sinSqr2 * sinSqr);
      const f32 sinSqr4 = (sinSqr3 * sinSqr);
      const f32 sinSqr5 = (sinSqr4 * sinSqr);

      const f32 p1 = (-0.045423f * sinSqr4) 
         + (0.74990f * sinSqr3) 
         - (3.24684f * sinSqr2) 
         + (5.46374f * sinSqr) 
         - 3.25565f;

      const f32 p2 = sinSqr5
         - (1.34515f * sinSqr4) 
         + (8.92015f * sinSqr3) 
         - (26.73322f * sinSqr2) 
         + (37.05174f * sinSqr) 
         - 19.53394f;

      return signal * (sin + (sin * sinSqr * (p1 * (HInvF(p2, 2)))));
}

f64 SLog2(f64 number)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;
      
      const f64 fExp = (f64) ((reInterpret.uN >> 52) & 0x7FF) - 1023.0;

      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFF) | 0x3FF0000000000000;

      const f64 x = (reInterpret.fN - 1.0) * HInv(reInterpret.fN + 1.0, 2);
      const f64 xSqr = x * x;
      
      const f64 mExp = x * (((0.412198511654412 * xSqr + 0.577078016299787) * xSqr + 0.961796693126021) * xSqr + 2.885390081777927);

      return fExp + mExp;
}

f32 SLog2F(f32 number)
{
      U32F32ToF32U32 reInterpret;
      reInterpret.fN = number;
      
      const f32 fExp = (f32) ((reInterpret.uN >> 23) & 0xFF) - 127.0f;

      reInterpret.uN = (reInterpret.uN & 0x7FFFFF) | 0x3F800000;
      
      const f32 x = reInterpret.fN - 1.0f; 

      const f32 mExp = (((-0.086566724f * x + 0.231248074f) * x - 0.334169527f) * x + 0.470030026f) * x - 0.718274794f;
      
      return (fExp + 1.0f) + (mExp * x);
}

f64 SLogN(f64 number)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;
      
      const f64 fExp = (f64) ((reInterpret.uN >> 52) & 0x7FF) - 1023.0;

      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFF) | 0x3FF0000000000000;

      const f64 x = (reInterpret.fN - 1.0) * HInv(reInterpret.fN + 1.0, 2);
      const f64 xSqr = x * x;
      
      const f64 mExp = x * (((0.412198511654412 * xSqr + 0.577078016299787) * xSqr + 0.961796693126021) * xSqr + 2.885390081777927);

      return (fExp + mExp) * 0.693147180559945309417;
}

f32 SLogNF(f32 number)
{
      U32F32ToF32U32 reInterpret;
      reInterpret.fN = number;
      
      const f32 fExp = (f32) ((reInterpret.uN >> 23) & 0xFF) - 127.0f;

      reInterpret.uN = (reInterpret.uN & 0x7FFFFF) | 0x3F800000;
      
      const f32 x = reInterpret.fN - 1.0f; 

      const f32 mExp = (((-0.086566724f * x + 0.231248074f) * x - 0.334169527f) * x + 0.470030026f) * x - 0.718274794f;
      
      return ((fExp + 1.0f) + (mExp * x)) * 0.69314718f;
}

f64 SLog10(f64 number)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;
      
      const f64 fExp = (f64) ((reInterpret.uN >> 52) & 0x7FF) - 1023.0;

      reInterpret.uN = (reInterpret.uN & 0x000FFFFFFFFFFFFF) | 0x3FF0000000000000;

      const f64 x = (reInterpret.fN - 1.0) * HInv(reInterpret.fN + 1.0, 2);
      const f64 xSqr = x * x;
      
      const f64 mExp = x * (((0.412198511654412 * xSqr + 0.577078016299787) * xSqr + 0.961796693126021) * xSqr + 2.885390081777927);

      return (fExp + mExp) * 0.301029995663981195214;
}

f32 SLog10F(f32 number)
{
      U32F32ToF32U32 reInterpret;
      reInterpret.fN = number;
      
      const f32 fExp = (f32) ((reInterpret.uN >> 23) & 0xFF) - 127.0f;

      reInterpret.uN = (reInterpret.uN & 0x7FFFFF) | 0x3F800000;
      
      const f32 x = reInterpret.fN - 1.0f; 

      const f32 mExp = (((-0.086566724f * x + 0.231248074f) * x - 0.334169527f) * x + 0.470030026f) * x - 0.718274794f;
      
      return ((fExp + 1.0f) + (mExp * x)) * 0.30103f;
}