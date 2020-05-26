/******************************************************************************
FFT²
*******************************************************************************/
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#define PI              3.1415927
#define ADC_DEF         1024
#define SMP_N       16
#define RE(n,b) *(n+(b<<1))
#define I(n) *(n+1)
#define IM(n,b) *(n+(b<<1)+1)
void
fftf_rdx2 (short *sig, float *TF, int N, float *twiddles)
{
  int i, j, tmp_N, k, btf, lg2 = 0;
  short *sig_t;
  float tmp_r, tmp_i;
  float *mid, *twid, *base;

  //DFT2

  for (i = N; i; i >>= 1)
    lg2++;

  // Toutes les DFT 2 sont effectuées avec des additions réelles
  for (i = 0; i < N; i += 2)
    {
      sig_t = sig + i;
      *(TF + (i << 1)) = *sig_t + *(sig_t + 1);
      *(TF + (i << 1) + 2) = *sig_t - *(sig_t + 1);
    }



  for (tmp_N = 4, i = 1; i < lg2; i++, tmp_N <<= 1)
    {

      for (j = 0; j < N; j += tmp_N)
	{
	  /*Adresse de X[0+j] 
	     et de X[N/2+j]
	   */
	  base = TF + (j << 1);
	  mid = base + tmp_N;

	  tmp_N >>= 1;

	  k = (SMP_N) / tmp_N;

	  tmp_r = *(mid);
	  tmp_i = I (mid);

	  I (mid) = I (base) - tmp_i;
	  *(mid) = *(base) - tmp_r;

	  I (base) = I (base) + tmp_i;
	  *(base) = *(base) + tmp_r;

	  for (btf = 1; btf < tmp_N; btf++)
	    {
	      twid = twiddles + k * btf;
	      // Multiplication du twiddle
	      tmp_r = RE (mid, btf) * *(twid) - IM (mid, btf) * I (twid);
	      tmp_i = RE (mid, btf) * I (twid) + IM (mid, btf) ** (twid);

	      RE (mid, btf) = RE (base, btf) - tmp_r;
	      IM (mid, btf) = IM (base, btf) - tmp_i;

	      RE (base, btf) = RE (base, btf) + tmp_r;
	      IM (base, btf) = IM (base, btf) + tmp_i;
	      twid += k;
	    }
	  tmp_N <<= 1;

	}
    }
}

void
fftf_rdx4 (short *sig, float *TF, int N, float *twiddles)
{
  int i, j, k, tmp_N, btf, lg4 = 0;
  short *sig_t, tmp;
  float bf[8], tmp_r, tmp_i;
  float *mid, *twid, *pos[4];


  for (i = N; i; i >>= 2)
    lg4++;

//DFT 4
  for (i = 0; i < N; i += 4)
    {
      *pos = TF + (i << 1);
      tmp = *sig + i;

      **pos = *(sig + i) + *(sig + i + 1) + *(sig + i + 2) + *(sig + i + 3);

      RE (*pos, 1) = *(sig + i) - *(sig + i + 2);
      IM (*pos, 1) = -*(sig + i + 3) + *(sig + i + 1);

      RE (*pos, 2) =
	*(sig + i) - *(sig + i + 1) + *(sig + i + 2) - *(sig + i + 3);

      RE (*pos, 3) = RE (*pos, 1);
      IM (*pos, 3) = -IM (*pos, 2);

    }
  //FFT(4^i+1) 
  for (tmp_N = 16, i = 1; i < 2; i++, tmp_N <<= 2)
    {
  k = (SMP_N / tmp_N)<<1;

      //Sur l'ensemble du signal
      for (j = 0; j < N; j += tmp_N)
	{
	  *pos = TF + (j << 1);

	  *bf = **pos;

	  // Premier butterfly à pour twiddle factor 1
	  for (i = 1; i < 4; i++)
	    {
	      pos[i] = pos[i - 1] + (tmp_N >> 1);
	      RE (bf, i) = *pos[i];
	      IM (bf, i) = I (pos[i]);
	    }

	  //Le 1er terme et le N/2 terme sont réels purs
	  **pos = *bf + RE (bf, 1) + RE (bf, 2) + RE (bf, 3);
	  *pos[2] = *bf - RE (bf, 1) + RE (bf, 2) - RE (bf, 3);

	  *pos[1] = *bf + IM (bf, 1) - RE (bf, 2) - IM (bf, 3);
	  I (pos[1]) = RE (bf, 3) - RE (bf, 1) - IM (bf, 2);

	  *pos[3] = *pos[1];
	  I (pos[3]) = -I (pos[1]);



	  for (btf = 1; btf < tmp_N >> 2; btf++)
	    {

	      *pos = TF + (j << 1) + (btf << 1);

	      for (i = 1; i < 4; i++)
		pos[i] = pos[i - 1] + (tmp_N >> 1);

	 	  printf("\n");

		  RE (bf, 0) = *pos[0] ;
		  IM (bf, 0) = I(pos[0]);
	      for (i = 1; i < tmp_N >> 2; i++)
		{
/*		  printf ("\n i * k * btf %d ", i * k * btf);
*/		  twid = twiddles + (i *btf  * k);

		 /* printf("%f ",*pos[i] );
		  printf ("\nTWEED: %f\n\t%f", *twid, I (twid));
*/
		  RE (bf, i) =*pos[i] * *(twid) - I(pos[i]) * I (twid);
		  IM (bf, i) = *pos[i] * I (twid) + I (pos[i]) ** (twid);
		}
		//multiplier par la matrice W4;
//		*pos[1] =  *(bf)-RE(bf,2)   ;
//		I(pos[1]) = I(bf) ;
		
		if(btf == 1){

		    	*pos[0]=*(bf)+RE(bf,3) ; 
			I(pos[0])=I(bf)+IM(bf,1)+IM(bf,2)+IM(bf,3) ; 
			for(i=0;i<4;i++)
			printf("\n%d : RE : %f Im : %f",i,RE(bf,i),IM(bf,i));
			
		}
		if(btf == 2){
		    	*pos[0]=*(bf)+RE(bf,1)+RE(bf,2)+RE(bf,3) ; 
			I(pos[0])=I(bf)+IM(bf,1)+IM(bf,2)+IM(bf,3) ; 

			*pos[2]=
                                 *bf -RE(bf,1)+RE(bf,2) - IM(bf,3);
			I(pos[2])= 
                                 I(bf) -IM(bf,1)+IM(bf,2) - IM(bf,3)  ;

			*pos[1] =*pos[2];		 
			I(pos[1]) = -I(pos[2]);
			
			*pos[3]=*pos[0];
			I(pos[3]) = -I(pos[0]);
		}
/*         switch(i & 0x3){
                      case 0 : 
                                tmp_r = *bf +RE(bf,1)+RE(bf,2) + RE(bf,3);
                                tmp_i = I(bf) +IM(bf,1)+IM(bf,2) + IM(bf,3)  ;
                      break;
                        case 1 : 
                                tmp_r = *bf +IM(bf,1)-RE(bf,2) + IM(bf,3);
                                tmp_i = I(bf) +RE(bf,1)-IM(bf,2) + RE(bf,3)  ;
                      break;
                        case 2 : 
                                tmp_r = *bf -RE(bf,1)+RE(bf,2) - IM(bf,3);
                                tmp_i = I(bf) -IM(bf,1)+IM(bf,2) - IM(bf,3)  ;
                      break;
                        case 3 : 
                                tmp_r = *bf   -IM(bf,1)   -RE(bf,2) + IM(bf,3);
                                tmp_i = I(bf) -RE(bf,1) -IM(bf,2) + IM(bf,3)  ;
                      break;
                        default:
                        break;
*/


	    }
	}
    }
}

     // Permet d'obtenir les twiddles factors
     // Pour une fft sur b termes
float *
get_twiddles_rdx4 (int b)
{

  int N = b << 1;

  float *twiddles = (float *) malloc (sizeof (float) * b);
  int half = N >> 1, ind;

  for (ind = 0; ind < b; ind += 2)
    {
      *(twiddles + ind) = cos (2*ind * PI/N );
      *(twiddles + ind + 1) = -sin (2*ind * PI /N );
    }

  return twiddles;
}
     // Permet d'obtenir les twiddles factors
     // Pour une fft sur b termes
float *
get_twiddles (int b)
{

  int N = b << 1;

  float *twiddles = (float *) malloc (sizeof (float) * b);
  int half = N >> 1, ind;

  for (ind = 0; ind < b; ind += 2)
    {
      *(twiddles + ind) = cos (ind * PI /half);
      *(twiddles + ind + 1) = -sin (ind * PI / half);
    }

  return twiddles;
}

     //Tableau de renversement
static unsigned char br_4[] = {
  0x0, 0x8, 0x4, 0xC, 0x2, 0xA, 0x6, 0xE,
  0x1, 0x9, 0x5, 0xD, 0x3, 0xB, 0x7, 0xF
};

     // Renverse les bits d'un short pour radix 2
int
rv_16_rdx2 (unsigned int j)
{
  return (0x0000ffff | (br_4[(j & 0x000f)] << 12) |
	  (br_4[(j & 0x00f0) >> 4] << 8) |
	  (br_4[(j & 0x0f00) >> 8] << 4) | (br_4[(j & 0x0f000) >> 12]));
}

     // Renverse un tableau de short pour radix 2
void
rvs_16_rdx2 (int *sig)
{
  unsigned int i;
  int tm[SMP_N];
  for (i = 0; i < SMP_N; i++)
    // i renversé
    tm[i] = sig[rv_16_rdx2 (i)];

  for (i = 0; i < SMP_N; i++)
    sig[i] = tm[i];

}


     // Renverse un tableau de short pour radix 4
void
rvs_16_rdx4 (short *sig, int N)
{
  int i, tmp_N = N >> 2;

  int tm[SMP_N];
  for (i = 0; i < tmp_N; i++)
    {
      tm[(i << 2)] = sig[i];
      tm[(i << 2) + 1] = sig[(i) + (tmp_N)];
      tm[(i << 2) + 2] = sig[(i) + (tmp_N << 1)];
      tm[(i << 2) + 3] = sig[(i) + (tmp_N << 1) + tmp_N];
    }

  for (i = 0; i < N; i++)
    sig[i] = tm[i];
}

int
main ()
{
  int i, N = SMP_N, b = N << 1;

  /* Tableau de complexes de taille N */
  float *TF = (float *) calloc (b, sizeof (float) * b);
  short *sig = malloc (sizeof (short) * (N));
  short test[16] =
    { 15, 7, 2, 9, 21, 15, 23, 10, 30, 14, 5, 7, 15, 14, 20, 16 };
  //float *twiddles = get_twiddles (N );
  float *twiddles = get_twiddles_rdx4(N);
  struct timespec now, bf;
  if (!sig || !twiddles || !TF)
    return 0;

  printf ("TF sur %d termes\n", N);
  printf ("\n\tBesoins mémoires :");
  printf ("\nTF :\t %d o \n", b * 4);
  printf ("SIG : \t%d o \n", N * 2);

  // Tablieau renversant le tableau du signal

  for (i = 0; i < N; i++)
    *(sig + i) = test[i];


//  rvs_16_rdx2 (sig);
  timespec_get (&bf, TIME_UTC);
//  fftf_rdx2 (sig, TF, N, twiddles);
  timespec_get (&now, TIME_UTC);
  printf ("\nRADIX 2\nTemps écoulé :\t %lf ms\n",
	  (now.tv_nsec - bf.tv_nsec) / 1e3);


  rvs_16_rdx4 (sig, N);

  for (i = 0; i < N; i++)
    printf ("\n%d = %d", i, sig[i]);

  timespec_get (&bf, TIME_UTC);
  fftf_rdx4 (sig, TF, N, twiddles);

  timespec_get (&now, TIME_UTC);
  printf ("\nRADIX 4\nTemps écoulé :\t %*lf ms\n", 4,
	  (now.tv_nsec - bf.tv_nsec) / 1e3);

  for (int i = 0; i < b; i += 2)
    {
      printf ("%d R %-2.2f\t\t", i / 2, *(TF + i));
      printf ("%-2.2f I\n", I (TF + i));
    }

  free (TF);
  free (twiddles);
  free (sig);
  return 0;
}
