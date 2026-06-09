#include "numbers.h"

#include <stdio.h>

void NormalizeFloatNumber(double number, double* pNormalized, i32* pExponent, i32* pSignal)
{
      const u64 numberBytes = GetFloatNumberBytes(number);
      int hasCorrection = 0;

      u64 numberMantissa = GetFloatNumberMantissa(numberBytes);
      i32 numberSignal = GetFloatNumberSignal(numberBytes);
      i32 numberExponent = GetFloatNumberExponent(numberBytes);

      if ((numberExponent % 2) != 0)
      {
            numberExponent -= 1;
            hasCorrection = 1;
      }

      double normalizedNumber = CreateFloatNumber(numberMantissa, 1023, 0);

      *pNormalized = normalizedNumber * (hasCorrection ? 2.0 : 1.0);
      *pExponent = numberExponent;
      *pSignal = numberSignal;
}