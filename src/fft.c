/******************************************************************************
			PROJET 	S6
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
*******************************************************************************/
#include <stdlib.h>
#include <math.h>
#define PI 3.1415927
#define X3(n) ((n<<1) + n)
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
	  k = (N) / tmp_N;
	  tmp_r = *(mid);
	  tmp_i = I (mid);

	  I (mid) = I (base) - tmp_i;
	  *(mid) = *(base) - tmp_r;

	  I (base) = I (base) + tmp_i;
	  *(base) = *(base) + tmp_r;

	  for (btf = 1; btf < tmp_N; btf++)
	    {
	      //W btf 
	      twid = twiddles + k * btf;

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
	//W0
	*TF =(float) *(sig ) + *(sig +  1) + *(sig +  2) + *(sig +  3);
	//W1
      	RE (TF, 1) = (float)*(sig ) - *(sig  + 2);
      	IM (TF, 1) = (float)*(sig  + 3) - *(sig  + 1);
	//W2
      	RE (TF, 2) =(float)
	*(sig ) - *(sig +  1) + *(sig  + 2) - *(sig  + 3);
	//W3
      	RE (TF, 3) = (float)RE (TF, 1);
      	IM (TF, 3) = (float)*(sig +  1) - *(sig + 3);
//FFT 4
  for (i = 4; i < N; i += 4)
    {
      *pos = TF + (i << 1);
//W0
      **pos = *(sig + i) + *(sig + i + 1) + *(sig + i + 2) + *(sig + i + 3);
//W1
      RE (*pos, 1) = *(sig + i) - *(sig + i + 2);
      IM (*pos, 1) = *(sig + i + 3) - *(sig + i + 1);
//W2
      RE (*pos, 2) =
	*(sig + i) - *(sig + i + 1) + *(sig + i + 2) - *(sig + i + 3);
//W3
      RE (*pos, 3) = RE (*pos, 1);
      IM (*pos, 3) = *(sig + i + 1) - *(sig + i + 3);
    }


//FFT(4^(lvl+1)) 
  for (tmp_N = 16, lvl = 1; lvl < lg4; lvl++, tmp_N <<= 2)
    {
      k = N / tmp_N;

//Sur l'ensemble du signal
      for (j = 0; j < N; j += tmp_N)
	{
	  //Position de FFT courante
	  *pos = TF+ (j<<1);
	 //  + (j *2);
	  //N/4 (TF*2 termes)
	  tmp_N >>= 1;

	  //Récupération N/4,N/2,3N/4
	  for (i = 1; i < 4; i++)
	    pos[i] = pos[i - 1] + (tmp_N);
	  //Calcul du premier papillon
	  for (i = 0; i < 4; i++){
	      RE (bf, i) = *pos[i];
	      IM(bf,i) = I(pos[i]);
	}

	      //X(0) 
	      **pos = W_0_4R (bf);
	      I (*pos) = W_0_4I (bf);
      //X(N/4)
	      *pos[1] = W_1_4R (bf);
	      I (pos[1]) = W_1_4I (bf);

	      //X(N/2)ne peut être obtenu par symétrie
	      *pos[2] = W_2_4R (bf);
	      I (pos[2]) = W_2_4I (bf);

	      for (btf = 1; btf < tmp_N >> 1; btf++)
		{
		    
	  //Haut du papillon
		  *pos = TF + (j*2) + (btf << 1);
		  *bf = *pos[0];
		  I (bf) = I (pos[0]);
		  //Récupération N/4,N/2,3N/4
		  for (i = 1; i < 4; i++)
		    pos[i] = pos[i - 1] + (tmp_N);

		  for (i = 1; i < 4; i++)
		    {
		      twid = twiddles + (i * (btf) * k * 2);
		      //W (i * btf) * N i/4
		      RE (bf, i) = *pos[i] * *(twid) - I (pos[i]) * I (twid);
		      IM (bf, i) = *pos[i] * I (twid) + I (pos[i]) ** (twid);
		    }
		  //X(4k) 
		  *pos[0] = W_0_4R (bf);
		  I (pos[0]) = W_0_4I (bf);
		  //X(4k+1) 
		  *pos[1] = W_1_4R (bf);
		  I (pos[1]) = W_1_4I (bf);
		}
	      //Les DFTS réelles sont symétriques
	      *pos = TF + (j* 2) + (tmp_N << 1);
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

// Pour une fft sur b termes
float *
get_twiddles (int b)
{
  int ind,max = b<<1;
  float *twiddles = (float *) malloc (sizeof (float) * (max));

  for (ind = 0; ind < max; ind += 2)
    {
      //e -i2pi/N 
      *(twiddles + ind) = cos (ind * PI / b);
      I (twiddles + ind) = -sin (ind * PI / b);
    }
  return twiddles;
}

