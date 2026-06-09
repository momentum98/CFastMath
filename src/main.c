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

      QueryPerformanceCounter(&start);

      //printf("%f\n", HardwareSqrt(12308));
      
      QueryPerformanceCounter(&end);
      
      double delay = (double) (end.QuadPart - start.QuadPart) * 1000.0 / freq.QuadPart;

      printf("Time: %.9f\n", delay);
}