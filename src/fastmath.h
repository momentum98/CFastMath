#ifndef FAST_MATH_H
#define FAST_MATH_H

typedef unsigned long long u64;
typedef unsigned int       u32;
typedef long long          i64;
typedef int                i32;

double SoftwareSqrtFrstLvl(double number, int iterations);
float SoftwareSqrtFrstLvlF(float number, int iterations);
double SoftwareSqrtScndLvl(double number, int iterations);
float SoftwareSqrtScndLvlF(float number, int iterations);

double HardwareSqrt(double number);
float HardwareSqrtF(float number);
float HardwareFastInverseSqrtF(float number);
float HardwareFastSqrtF(float number);

#endif