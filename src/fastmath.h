#ifndef FAST_MATH_H
#define FAST_MATH_H

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

f64 WrapTo180(f64 angle);
f32 WrapTo180F(f32 angle);

f64 WrapTo360(f64 angle);
f32 WrapTo360F(f32 angle);

f64 RadToDeg(f64 rad);
f32 RadToDegF(f32 rad);

f64 RadToTurn(f64 rad);
f32 RadToTurnF(f32 rad);

f64 DegToRad(f64 deg);
f32 DegToRadF(f32 deg);

f64 DegToTurn(f64 deg);
f32 DegToTurnF(f32 deg);

f64 TurnToRad(f64 turn);
f32 TurnToRadF(f32 turn);

f64 TurnToDeg(f64 turn);
f32 TurnToDegF(f32 turn);

f64 SSqrt(f64 number, i32 iterations);
f32 SSqrtF(f32 number, i32 iterations);

f64 SSqrt2(f64 number, i32 iterations);
f32 SSqrt2F(f32 number, i32 iterations);

f64 HSqrt(f64 number);
f32 HSqrtF(f32 number);

f32 HFInvSqrtF(f32 number);

f32 HFSqrtF(f32 number);

f64 SExp(f64 number, f64 exp);
f32 SExpF(f32 number, f32 exp);

f64 SFExp(f64 number, f64 exp);
f32 SFExpF(f32 number, f32 exp);

f64 SMax(f64 a, f64 b);
f32 SMaxF(f32 a, f32 b);

f64 HMax(f64 a, f64 b);
f32 HMaxF(f32 a, f32 b);

f64 SMin(f64 a, f64 b);
f32 SMinF(f32 a, f32 b);

f64 HMin(f64 a, f64 b);
f32 HMinF(f32 a, f32 b);

f64 SClamp(f64 value, f64 min, f64 max);
f32 SClampF(f32 value, f32 min, f32 max);

f64 HClamp(f64 value, f64 min, f64 max);
f32 HClampF(f32 value, f32 min, f32 max);

f64 HRound(f64 number);
f32 HRoundF(f32 number);

f64 HFloor(f64 number);
f32 HFloorF(f32 number);

f64 HCeil(f64 number);
f32 HCeilF(f32 number);

u32 IsNegative(f64 number);
u32 IsNegativeF(f32 number);

u32 SSign(f64 number);
u32 SSignF(f32 number);

f64 SAbs(f64 number);
f32 SAbsF(f32 number);

f64 SCopySign(f64 number, f64 signNumber);
f32 SCopySignF(f32 number, f32 signNumber);

f64 SFLerp(f64 a, f64 b, f64 p);
f32 SFLerpF(f32 a, f32 b, f32 p);

f64 SLerp(f64 a, f64 b, f64 p);
f32 SLerpF(f32 a, f32 b, f32 p);

f64 SEaseIn(f64 a, f64 b, f64 p);
f32 SEaseInF(f32 a, f32 b, f32 p);

f64 SEaseOut(f64 a, f64 b, f64 p);
f32 SEaseOutF(f32 a, f32 b, f32 p);

f64 SEaseInOut(f64 a, f64 b, f64 p);
f32 SEaseInOutF(f32 a, f32 b, f32 p);

f64 SInv(f64 number, u32 iterations);
f32 SInvF(f32 number, u32 iterations);

f64 HInv(f64 number, u32 iterations);
f32 HInvF(f32 number);

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
f64 SLog2F(f64 number);

f64 SLogN(f64 number);
f64 SLogNF(f64 number);

f64 SLog10(f64 number);
f64 SLog10F(f64 number);

#endif