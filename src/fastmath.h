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

f64 HRound(f64 number);
f32 HRoundF(f32 number);

f64 HFloor(f64 number);
f32 HFloorF(f32 number);

f64 HCeil(f64 number);
f32 HCeilF(f32 number);

#endif