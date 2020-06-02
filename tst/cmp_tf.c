/******************************************************************************
			PROJET 	S6
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
*******************************************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

float
cmp_tf (float *ref, float *tf2, int N)
{
  int i;
  float f = 0;
  for (i = 0; i < N * 2; i++)
    {
      if (*(ref + i) != *(tf2 + i))
	f += fabs (*(ref + i) / *(tf2 + i) - 1);

    }
  return f /= N;

}
