#include "fastmath.h"

#include <stdio.h>
#include <windows.h>
#include <math.h>

int main(int argc, char** args)
{
      LARGE_INTEGER start;
      LARGE_INTEGER end;
      LARGE_INTEGER freq;

      QueryPerformanceFrequency(&freq);

      /*
      printf("%f\n", SqrtSecondLvl(12308, 1));
      */

      double delay = (double) (start.QuadPart - end.QuadPart) * 1000.0 / freq.QuadPart;
      
      QueryPerformanceCounter(&start);
      
      printf("Time: %.9f", delay);
}