#include "fastmath.h"

#include "numbers.h"

#include <stdlib.h>
#include <stdio.h>

const static double sqrtLUTStepSize = 0.1f;
const static double inverseSqrtLUTStepSize = 1 / sqrtLUTStepSize;

const static double sqrtLUT[41] = 
{
      0.0,
      0.3162277660168379,
      0.4472135954999579,
      0.5477225575051661,
      0.6324555320336758,
      0.7071067811865476,
      0.7745966692414834,
      0.8366600265340756,
      0.8944271909999159,
      0.9486832980505138,
      1.0,
      1.0488088481701516,
      1.0954451150103323,
      1.140175425099138,
      1.1832159566199232,
      1.224744871391589,
      1.2649110640673518,
      1.3038404810405297,
      1.3416407864998738,
      1.378404875209022,
      1.4142135623730951,
      1.449137674618944,
      1.4832396974191326,
      1.51657508881031,
      1.5491933384829668,
      1.5811388300841898,
      1.61245154965971,
      1.6431676725154984,
      1.6733200530681511,
      1.7029386365926402,
      1.7320508075688772,
      1.7606816861038855,
      1.7888543819998317,
      1.816590212458495,
      1.8439088914585775,
      1.8708286933869707,
      1.8973665961010276,
      1.9235384061671343,
      1.9493588689617927,
      1.9748417658131499,
      2.0
};

static void GetClosestSqrtLUTData(double number, double* pLUTValue, double* pLUTSqrValue)
{
      i32 index = (int)((number * inverseSqrtLUTStepSize) + 0.5);

      *pLUTValue    = sqrtLUT[index];
      *pLUTSqrValue = sqrtLUTStepSize * index;
}

static double SqrtIterateFirstLvl(double number, double try, double trySqr)
{
      const double correction = (number - trySqr) / (2.0 * try);
      return try + correction;
}

double SqrtFirstLvl(double number, int iterations)
{ 
      if (number < 0)
      {
            return -1;
      }
      
      if (number == 0 || number == 1)
      {
            return number;
      }

      i32 signal = -1;
      i32 exponent = -1;

      double normalizedNumber = -1.0;

      NormalizeFloatNumber(number, &normalizedNumber, &exponent, &signal);

      double LUTValue = -1.0;
      double LUTSqrValue = -1.0;

      GetClosestSqrtLUTData(normalizedNumber, &LUTValue, &LUTSqrValue);

      double try = LUTValue;
      double trySqr = LUTSqrValue;

      for (i32 i = 0; i < iterations; i++)
      {
            try = SqrtIterateFirstLvl(normalizedNumber, try, trySqr);
            trySqr = try * try;
      }

      u64 tryBytes = GetFloatNumberBytes(try);
      u64 tryMantissa = GetFloatNumberMantissa(tryBytes);

      return CreateFloatNumber(tryMantissa, (u32) ((exponent >> 1) + 1023), signal);
}

static double SqrtIterateSecondLvl(double number, double try, double trySqr)
{
      return try * ((try * try) + (3 * number)) / ((3 * (try * try)) + number);
}

double SqrtSecondLvl(double number, int iterations)
{ 
      if (number < 0)
      {
            return -1;
      }
      
      if (number == 0 || number == 1)
      {
            return number;
      }

      i32 signal = -1;
      i32 exponent = -1;

      double normalizedNumber = -1.0;

      NormalizeFloatNumber(number, &normalizedNumber, &exponent, &signal);

      double LUTValue = -1.0;
      double LUTSqrValue = -1.0;

      GetClosestSqrtLUTData(normalizedNumber, &LUTValue, &LUTSqrValue);

      double try = LUTValue;
      double trySqr = LUTSqrValue;

      for (i32 i = 0; i < iterations; i++)
      {
            try = SqrtIterateSecondLvl(normalizedNumber, try, trySqr);
            trySqr = try * try;
      }

      u64 tryBytes = GetFloatNumberBytes(try);
      u64 tryMantissa = GetFloatNumberMantissa(tryBytes);

      return CreateFloatNumber(tryMantissa, (u32) ((exponent >> 1) + 1023), signal);
}