#ifndef FAST_MATH_H
#define FAST_MATH_H

typedef unsigned long long u64;
typedef unsigned int       u32;
typedef long long          i64;
typedef int                i32;
typedef float              f32;
typedef double             f64;

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

f64 HRound(f64 number);
f32 HRoundF(f32 number);

f64 HFloor(f64 number);
f32 HFloorF(f32 number);

f64 HCeil(f64 number);
f32 HCeilF(f32 number);

f64 SAbs(f64 number);
f32 SAbsF(f32 number);

f64 SCopySign(f64 number, f64 signNumber);
f32 SCopySignF(f32 number, f32 signNumber);

f64 SFLerp(f64 a, f64 b, f64 p);
f32 SFLerpF(f32 a, f32 b, f32 p);

f64 SLerp(f64 a, f64 b, f64 p);
f32 SLerpF(f32 a, f32 b, f32 p);

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

#endif