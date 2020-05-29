
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

  int i, ind, n_bit = -1;
  short tmp;

  for (i = N; i; i >>= 1)
    n_bit++;

  ind = (br_4[(i & 0x000f)] << 12);
//
  for (i = 1; i < N - 1; i++)
    {
      tmp = sig[i];

      ind = i;
      sig[i] = sig[ind];
      sig[ind] = tmp;

// i renversé
    }
}

// Renverse un tableau de short à indice renversé pour radix 4
void
rvs_16_rdx4 (short *sig, int N)
{

  int  i  ;

int j = 0, N1,N2 =N<<1 ;

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
      while (j >= 3*N1)
	{
	  j -= 3*N1;
	  N1 >>= 2;
	}
      j += N1;
    }
}
