/******************************************************************************
			PROJET 	S6
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
*******************************************************************************/
#define X3(n) ((n<<1) + n)

//Tableau de renversement
static unsigned char br_4[] = {
  0x0, 0x8, 0x4, 0xC, 0x2, 0xA, 0x6, 0xE,
  0x1, 0x9, 0x5, 0xD, 0x3, 0xB, 0x7, 0xF
};


// Renverse un tableau de short a indice renversé pour radix 2

void
rvs_16_rdx2 (short *sig, int N)
{

  int i;
  unsigned int tmp, r, s, lg2 = 0;
  short tm[N];

//
  switch (N)
    {
    case 1 << 4:
      for (i = 0; i < N; i++)
	{
	  tm[i] = sig[br_4[i]];
	}
      for (i = 0; i < N; i++)
	{
	  sig[i] = tm[i];
	}

      break;
    case 1 << 8:

      for (i = 0; i < N; i++)
	{
	  tm[i] = sig[((br_4[(i & 0x0f)] << 4) | (br_4[(i & 0xf0) >> 4]))];
	}

      break;
    case 1 << 12:
      for (i = 0; i < N; i++)
	{
	  tm[i] = sig[br_4[(i & 0x0f)] << 8 |
		      br_4[(i & 0xf0) >> 4] << 4 | br_4[(i & 0xf00) >> 8]];
	}
      break;
    case 1 << 16:
      for (i = 0; i < N; i++)
	{
	  tm[i] = sig[br_4[(i & 0x0f)] << 12 |
		      br_4[(i & 0xf0) >> 4] << 8 |
		      br_4[(i & 0xf00) >> 8] << 4 | br_4[(i & 0xf000) >> 12]];
	}
      break;
    default:
      while (N >>= 1)
	lg2++;
      N = 2;
      N <<= (lg2 - 1);
      for (i = 1; i < N; i++)
	{

	  r = i;
	  s = (sizeof (i) * 8 - 1);
	  for (tmp = i >> 1; tmp; tmp >>= 1)
	    {
	      r <<= 1;
	      r |= tmp & 1;
	      s--;
	    }
	  r <<= s;
	  r >>= (sizeof (int) * 8) - lg2;
	  tm[i] = sig[r];
	}
      break;
    }


  for (i = 1; i < N - 1; i++)
    {
      sig[i] = tm[i];
    }
}

// Renverse un tableau de short à indice renversé pour radix 4
void
rvs_16_rdx4 (short *sig, int N)
{

  int i;

  int j = 0, N1;

  short tmp;
  for (i = 0; i < N - 1; i++)
    {
      if (i < j)
	{
	  tmp = sig[i];
	  sig[i] = sig[j];
	  sig[j] = tmp;
	}
      N1 = N >> 2;
      while (j >= 3 * N1)
	{
	  j -= 3 * N1;
	  N1 >>= 2;
	}
      j += N1;
    }
}
