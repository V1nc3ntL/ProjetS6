/******************************************************************************


FFT : radix 2 et 4


*******************************************************************************/


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>


#define PI 3.1415927


#define ADC_DEF 1024


#define SMP_N 16



#define RE(n,b) *(n+(b<<1))


#define I(n) *(n+1)


#define IM(n,b) *(n+(b<<1)+1)



#define W_0_4R(n) (RE(n,0)+RE(n,1)+RE(n,2)+RE(n,3))


#define W_0_4I(n) (IM(n,0)+IM(n,1)+IM(n,2)+IM(n,3))



#define W_1_4R(n) (*n+IM(n,1))-RE(n,2)-IM(n,3)


#define W_1_4I(n) (I(n)-RE(n,1)-IM(n,2)+RE(n,3))



#define W_2_4R(n) (RE(n,0)-RE(n,1)+RE(n,2)-RE(n,3))


#define W_2_4I(n) (IM(n,0)-IM(n,1)+IM(n,2)-IM(n,3))



#define W_3_4R(n) (RE(n,0)-IM(n,1)-RE(n,2)+IM(n,3))


#define W_3_4I(n) (IM(n,0)+RE(n,1)-IM(n,2)-RE(n,3))



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


  int lvl, i, j, k, tmp_N, btf, lg4 = 0;

  float bf[8];


  float *twid, *pos[4];



  for (i = N; i; i >>= 2)


    lg4++;


  lg4--;


//FFT 4
  for (i = 0; i < N; i += 4)


    {
      *pos = TF + (i << 1);

      **pos = *(sig + i) + *(sig + i + 1) + *(sig + i + 2) + *(sig + i + 3);

      RE (*pos, 2) =
	*(sig + i) - *(sig + i + 1) + *(sig + i + 2) - *(sig + i + 3);

      RE (*pos, 1) = *(sig + i) - *(sig + i + 2);
      IM (*pos, 1) = *(sig + i + 3) - *(sig + i + 1);

      RE (*pos, 3) = RE (*pos, 1) ;
      IM (*pos, 3) = *(sig + i + 1) - *(sig + i + 3);

    }


//FFT(4^(lvl+1)) 

  for (tmp_N = 16, lvl = 1; lvl < lg4; lvl++, tmp_N <<= 2)


    {

      k = SMP_N / tmp_N;

//Sur l'ensemble du signal

      for (j = 0; j < N; j += tmp_N)
	{
//Permet d'obtenir le quart de la TF
	  tmp_N >>= 1;

	  *pos = TF + (j << 1);

	  for (i = 1; i < 4; i++)
    pos[i] = pos[i - 1] + (tmp_N);

	  for (i = 0; i < 4; i++){
				  RE (bf, i) = *pos[i] ;

				  IM (bf, i) = I(pos[i]) ;
	  	}
	
//Calcul du premier papillon
      	  *pos[0] = W_0_4R (bf);
	  I (pos[0]) = W_0_4I (bf);

	  *pos[1] = W_1_4R (bf);
	  I (pos[1]) = W_1_4I (bf);

	  *pos[2] = W_2_4R (bf);
  I (pos[2]) = W_2_4I (bf);

	  *pos[3] = W_3_4R (bf);
	  I (pos[3]) = W_3_4I (bf);
	  

	  for (btf = 1; btf < tmp_N>>1; btf++)
	    {

	      *pos = TF + (j << 1) + (btf << 1);
	      *bf = *pos[0];
	      I (bf) = I (pos[0]);

	      for (i = 1; i < 4; i++)
		pos[i] = pos[i - 1] + (tmp_N);
		
	      	for (i = 1; i < 4; i++)
		{
		  twid = twiddles + (i * (btf) * k * 2);
		  RE (bf, i) = *pos[i] * *(twid) - I (pos[i]) * I (twid);
		  IM (bf, i) = *pos[i] * I (twid) + I (pos[i]) ** (twid);
		}
		
	     *pos[0] = W_0_4R (bf);
	     I(pos[0]) = W_0_4I (bf);

	      *pos[1] = W_1_4R (bf);
	      I (pos[1]) = W_1_4I (bf);

	    }
//Les DFTS réelles sont symétriques
	  *pos = TF + (j << 1) + (tmp_N << 1);
	  for (i = 1; i < tmp_N; i++)
	    {
	      RE (*pos, i) = RE (*pos, -i);
	      IM (*pos, i) = -IM (*pos, -i);
	    }
	// N courant
	  tmp_N <<= 1;
	}
    }
}



// Permet d'obtenir les twiddles factors


// Pour une fft sur b termes


float *
get_twiddles_rdx4 (int b)
{



  int N = b << 1;



  float *twiddles = (float *) calloc (b, sizeof (float) * b);


  int ind;



  for (ind = 0; ind < N; ind += 2)


    {


      *(twiddles + ind) = cos (ind * PI / b);


      I (twiddles + ind) = -sin (ind * PI / b);


    }



// *(twiddles + b ) = -1;



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


      *(twiddles + ind) = cos (ind * PI / half);


      I (twiddles + ind) = -sin (ind * PI / half);


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
/*
	tm[i] =sig[i<<2] ;
	tm[i+tmp_N] =sig[(i<<2)+1] ;
	tm[i+(tmp_N<<1)] =sig[(i<<2)+2] ;
	tm[i+(tmp_N<<1)+tmp_N] =sig[(i<<2)+3] ;
*/

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


  short test[64] = {
       	6, 13, 10, 16, 25, 14, 7, 3,
    	37,25,63,42,13,7,5,11,
   	5, 9, 17, 1, 2, 4, 12, 8,
    	7,21,17,12,11,10,4,21,
    	37,25,63,42,13,7,5,11,
   	5, 9, 17, 1, 2, 4, 12, 8,
       	6, 13, 10, 16, 25, 14, 7, 3,
    	7,21,17,12,11,10,4,21,


  };


//float *twiddles = get_twiddles (N );


  float *twiddles = get_twiddles_rdx4 (N);


  struct timespec now, bf;


  if (!sig || !twiddles || !TF)


    return 0;



  printf ("TF sur %d termes\n", N);


  printf ("\n\tBesoins mémoires :");


  printf ("\nTF :\t %d o \n", b * 4);


  printf ("SIG : \t%d o \n", N * 2);



  for (i = 0; i < N; i++)


    *(sig + i) =test[ i];



// rvs_16_rdx2 (sig);


  timespec_get (&bf, TIME_UTC);


// fftf_rdx2 (sig, TF, N, twiddles);


  timespec_get (&now, TIME_UTC);


  printf ("\nRADIX 2\nTemps écoulé :\t %lf ms\n",
	  (now.tv_nsec - bf.tv_nsec) / 1e3);



  rvs_16_rdx4 (sig, N);



  for (i = 0; i < N; i++)


    printf ("%d,", sig[i]);



  timespec_get (&bf, TIME_UTC);


  fftf_rdx4 (sig, TF, N, twiddles);



  timespec_get (&now, TIME_UTC);


  printf ("\nRADIX 4\nTemps écoulé :\t %*lf ms\n", 4,
	  (now.tv_nsec - bf.tv_nsec) / 1e3);



  for (int i = 0; i < b; i += 2)


    {


 //if(!(i%8)){

      printf ("%d R %+-2.2f\t\t", i / 2, *(TF + i));


      printf ("%+-2.2f I\n", I (TF + i));

 //}


    }



  free (TF);


  free (twiddles);


  free (sig);


  return 0;


}
