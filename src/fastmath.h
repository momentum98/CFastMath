#ifndef FAST_MATH_H
#define FAST_MATH_H

#include <immintrin.h>

typedef unsigned long long u64;
typedef unsigned int       u32;
typedef long long          i64;
typedef int                i32;
typedef float              f32;
typedef double             f64;

typedef struct
{
      f64 x;
      f64 y;
} Vector2;

typedef struct
{
      f32 x;
      f32 y;
} Vector2F;

typedef struct
{
      f64 x;
      f64 y;
      f64 z;
} Vector3;

typedef struct
{
      f32 x;
      f32 y;
      f32 z;
} Vector3F;

static inline f64 WrapTo180(f64 angle)
{
      return angle - 360.0 * (HFloor((angle + 180.0) * INV_FULL_ANGLE));
}

static inline f32 WrapTo180F(f32 angle)
{
      return angle - (360.0f * HFloorF((angle + 180.0f) * INV_FULL_ANGLEF));
}

static inline f64 WrapTo360(f64 angle)
{
      return angle - (360.0 * HFloor(angle * INV_FULL_ANGLE));
}

static inline f32 WrapTo360F(f32 angle)
{
      return angle - (360.0f * HFloorF(angle * INV_FULL_ANGLEF));
}

static inline f64 RadToDeg(f64 rad)
{
      return (rad * INV_PI) * 180.0;
}

static inline f32 RadToDegF(f32 rad)
{
      return (rad * INV_PIF) * 180.0f;
}

static inline f64 RadToTurn(f64 rad)
{
      return (rad * INV_PI * 0.5);
}

static inline f32 RadToTurnF(f32 rad)
{
      return (rad * INV_PIF * 0.5f);
}

static inline f64 DegToRad(f64 deg)
{
      return (deg * INV_FULL_ANGLE) * (2.0 * PI);
}

static inline f32 DegToRadF(f32 deg)
{
      return (deg * INV_FULL_ANGLEF) * (2.0f * PIF);
}

static inline f64 DegToTurn(f64 deg)
{
      return (deg * INV_FULL_ANGLE);
}

static inline f32 DegToTurnF(f32 deg)
{
      return (deg * INV_FULL_ANGLEF);
}

static inline f64 TurnToRad(f64 turn)
{
      return turn * (2.0 * PI);
}

static inline f32 TurnToRadF(f32 turn)
{
      return turn * (2.0f * PIF);
}

static inline f64 TurnToDeg(f64 turn)
{
      return turn * 360.0;
}

static inline f32 TurnToDegF(f32 turn)
{
      return turn * 360.0f;
}

static inline f64 HRound(f64 number)
{
      const __m128d reg = _mm_set_sd(number);
      const __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

static inline f32 HRoundF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}

static inline f64 HFloor(f64 number)
{
      const __m128d reg = _mm_set_sd(number);
      const __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

static inline f32 HFloorF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}

static inline f64 HCeil(f64 number)
{
      const __m128d reg = _mm_set_sd(number);
      const __m128d value = _mm_round_sd(reg, reg, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtsd_f64(value);
}

static inline f32 HCeilF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_round_ss(reg, reg, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

      return _mm_cvtss_f32(value);
}

static inline u32 IsNegative(f64 number)
{
      U64F64ToF64U64 reInterpret;
      reInterpret.fN = number;
      
      return (reInterpret.uN >> 63) & 1;
}

static inline u32 IsNegativeF(f32 number)
{
      U32F32ToF32U32 reInterpret;
      reInterpret.fN = number;
      
      return (reInterpret.uN >> 31) & 1;
}

static inline i32 SSign(f64 number)
{
      return (number > 0.0) - (number < 0.0);
}

static inline i32 SSignF(f32 number)
{
      return (number > 0.0f) - (number < 0.0f);
}

static inline f64 SAbs(f64 number)
{
      U64F64ToF64U64 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = (reInterpret.uN & 0x7FFFFFFFFFFFFFFF);

      return reInterpret.fN;
}

static inline f32 SAbsF(f32 number)
{
      U32F32ToF32U32 reInterpret;

      reInterpret.fN = number;
      reInterpret.uN = (reInterpret.uN & 0x7FFFFFFF);

      return reInterpret.fN;
}

static inline f64 SCopySign(f64 number, f64 signNumber)
{
      U64F64ToF64U64 reInterpretNum;
      U64F64ToF64U64 reInterpretSign;

      reInterpretSign.fN = signNumber;
      reInterpretNum.fN = number;
      
      reInterpretNum.uN = (reInterpretNum.uN & 0x7FFFFFFFFFFFFFFF) | (reInterpretSign.uN & 0x8000000000000000);

      return reInterpretNum.fN;
}

static inline f32 SCopySignF(f32 number, f32 signNumber)
{
      U32F32ToF32U32 reInterpretNum;
      U32F32ToF32U32 reInterpretSign;

      reInterpretSign.fN = signNumber;
      reInterpretNum.fN = number;

      reInterpretNum.uN = (reInterpretNum.uN & 0x7FFFFFFF) | (reInterpretSign.uN & 0x80000000);

      return reInterpretNum.fN;
}

static inline f64 SFLerp(f64 a, f64 b, f64 p)
{
      return a + (b - a) * p;
}

static inline f32 SFLerpF(f32 a, f32 b, f32 p)
{
      return a + (b - a) * p;
}

static inline f64 SLerp(f64 a, f64 b, f64 p)
{
      return (1 - p) * a + b * p;
}

static inline f32 SLerpF(f32 a, f32 b, f32 p)
{
      return (1 - p) * a + b * p;
}

static inline f64 SEaseIn(f64 a, f64 b, f64 p)
{
      const f64 finalProgress = (b - a) * (p * p);
      
      return a + finalProgress;
}

static inline f32 SEaseInF(f32 a, f32 b, f32 p)
{
      const f32 finalProgress = (b - a) * (p * p);
      
      return a + finalProgress;
}

static inline f64 SEaseOut(f64 a, f64 b, f64 p)
{
      const f64 finalProgress = (b - a) * (p * (2.0 - p));
      
      return a + finalProgress;
}

static inline f32 SEaseOutF(f32 a, f32 b, f32 p)
{
      const f32 finalProgress = (b - a) * (p * (2.0f - p));
      
      return a + finalProgress;
}

static inline f64 SEaseInOut(f64 a, f64 b, f64 p)
{
      const f64 finalProgress = (b - a) * (p * p * (3.0 - 2.0 * p));
      
      return a + finalProgress;
}

static inline f32 SEaseInOutF(f32 a, f32 b, f32 p)
{
      const f32 finalProgress = (b - a) * (p * p * (3.0f - 2.0f * p));
      
      return a + finalProgress;
}

f64 SInv(f64 number, u32 iterations);
f32 SInvF(f32 number, u32 iterations);

f64 HInv(f64 number, u32 iterations);
f32 HInvF(f32 number, u32 iterations);

f64 SSqrt(f64 number, u32 iterations);
f32 SSqrtF(f32 number, u32 iterations);

f64 SSqrt2(f64 number, u32 iterations);
f32 SSqrt2F(f32 number, u32 iterations);

static inline f64 HSqrt(f64 number)
{
      const __m128d reg = _mm_set_sd(number);
      const __m128d value = _mm_sqrt_sd(reg, reg);

      return _mm_cvtsd_f64(value);
}

static inline f32 HSqrtF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_sqrt_ss(reg);

      return _mm_cvtss_f32(value);
}

static inline f32 HFInvSqrtF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_rsqrt_ss(reg);

      return _mm_cvtss_f32(value);
}

static inline f32 HFSqrtF(f32 number)
{
      const __m128 reg = _mm_set_ss(number);
      const __m128 value = _mm_rsqrt_ss(reg);

      return _mm_cvtss_f32(value) * number;
}

f64 SExp(f64 number, f64 exp);
f32 SExpF(f32 number, f32 exp);

static inline f64 SMax(f64 a, f64 b)
{
      return ((a + b) + SAbs(a - b)) * 0.5;
}

static inline f32 SMaxF(f32 a, f32 b)
{
      return ((a + b) + SAbsF(a - b)) * 0.5f;
}

static inline f64 HMax(f64 a, f64 b)
{
      const __m128d aReg = _mm_set_sd(a);
      const __m128d bReg = _mm_set_sd(b);

      const __m128d value = _mm_max_sd(aReg, bReg);

      return _mm_cvtsd_f64(value);
}

static inline f32 HMaxF(f32 a, f32 b)
{
      const __m128 aReg = _mm_set_ss(a);
      const __m128 bReg = _mm_set_ss(b);

      const __m128 value = _mm_max_ss(aReg, bReg);

      return _mm_cvtss_f32(value);
}

static inline f64 SMin(f64 a, f64 b)
{
      return ((a + b) - SAbs(a - b)) * 0.5;
}

static inline f32 SMinF(f32 a, f32 b)
{
      return ((a + b) - SAbsF(a - b)) * 0.5f;
}

static inline f64 HMin(f64 a, f64 b)
{
      const __m128d aReg = _mm_set_sd(a);
      const __m128d bReg = _mm_set_sd(b);

      const __m128d value = _mm_min_sd(aReg, bReg);

      return _mm_cvtsd_f64(value);
}

static inline f32 HMinF(f32 a, f32 b)
{
      const __m128 aReg = _mm_set_ss(a);
      const __m128 bReg = _mm_set_ss(b);

      const __m128 value = _mm_min_ss(aReg, bReg);

      return _mm_cvtss_f32(value);
}

static inline f64 SClamp(f64 value, f64 min, f64 max)
{
      const f64 minPatched = SMax(value, min);
      const f64 maxPatched = SMin(minPatched, max);

      return maxPatched;
}

static inline f32 SClampF(f32 value, f32 min, f32 max)
{
      const f32 minPatched = SMaxF(value, min);
      const f32 maxPatched = SMinF(minPatched, max);

      return maxPatched;
}

static inline f64 HClamp(f64 value, f64 min, f64 max)
{
      const f64 minPatched = HMax(value, min);
      const f64 maxPatched = HMin(minPatched, max);

      return maxPatched;
}

static inline f32 HClampF(f32 value, f32 min, f32 max)
{
      const f32 minPatched = HMaxF(value, min);
      const f32 maxPatched = HMinF(minPatched, max);

      return maxPatched;
}

f64 SFExp(f64 number, f64 exp);
f32 SFExpF(f32 number, f32 exp);

f64 SCos(f64 turn);
f32 SCosF(f32 turn);

f64 SSin(f64 turn);
f32 SSinF(f32 turn);

f64 SFCos(f64 turn);
f32 SFCosF(f32 turn);

f64 SFSin(f64 turn);
f32 SFSinF(f32 turn);

f64 SFTan(f64 turn);
f32 SFTanF(f32 turn);

f64 STan(f64 turn);
f32 STanF(f32 turn);

f64 SAtan2(f64 x, f64 y);
f32 SAtan2F(f32 x, f32 y);

f64 SACos(f64 cos);
f32 SACosF(f32 cos);

f64 SASin(f64 sin);
f32 SASinF(f32 sin);

f64 SLog2(f64 number);
f32 SLog2F(f32 number);

f64 SLogN(f64 number);
f32 SLogNF(f32 number);

f64 SLog10(f64 number);
f32 SLog10F(f32 number);

#endif