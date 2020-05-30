/******************************************************************************
			PROJET 	S6
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
*******************************************************************************/


#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <stdlib.h>
#include "br_tst.h"
#include "fft_tst.h"
#include "string.h"
#define MAX 1<<14
int
main (int argc,char** argv)
{


  int i, N=4 ,tmp=0  ;
  char state=0;
  unsigned int * cmp_mul= (unsigned int *) calloc(1,sizeof(unsigned int)),* cmp_add=calloc(1,sizeof(unsigned int));
  struct timespec now, bf;
  /* Tableau de complexes de taille N */
  float *TF ;
  /* Tableau du signal */
  short *sig ;
  /* Tableau des twiddles */
  float *twiddles;


	if(argc==1){
 		for(N =16;N<=MAX;N<<=1){
			twiddles = get_twiddles (N);
			sig = calloc (N, sizeof (short) * (N));
			TF =   (float *) calloc (N<<1, sizeof (float) * (N<<1));
  			if (!sig || !twiddles || !TF|| !cmp_add || !cmp_mul)
    			{
      				printf ("Pas assez de mémoire\n");
			      return 0;
    			}

			printf ("\n\tTF sur %d termes ",N);

  			printf ("\nBesoins mémoires :");
 			printf ("\nTF :\t%d o \n", N * 8);
			printf ("TWI : \t%d o \n", N * 8);
			printf ("SIG : \t%d o \n", N * 2);
  			for (i = 0; i < N; i++)
			    *(sig + i) = i;

  			timespec_get (&bf, TIME_UTC);
  			rvs_16_rdx2 (sig,N);
			timespec_get (&now, TIME_UTC);
			printf ("\nRenversement de bits 2\nTemps écoulé :\t %lf ms\n",
			(now.tv_nsec - bf.tv_nsec) / 1e3);

  			timespec_get (&bf, TIME_UTC);
			fftf_rdx2 (sig, TF, N, twiddles,cmp_add,cmp_mul);
			timespec_get (&now, TIME_UTC);
			printf ("\n\tRADIX 2\nTemps écoulé :\t %lf ms",
			(now.tv_nsec - bf.tv_nsec) / 1e3);
			printf("\n%d additions complexes\n%d multiplications complexes\n",*cmp_add,*cmp_mul );
			    if(!state){
			*cmp_mul = 0;
			*cmp_add= 0 ;
  			for (i = 0; i < N; i++)
			    *(sig + i) = i;

  			timespec_get (&bf, TIME_UTC);
			rvs_16_rdx4 (sig, N);
			timespec_get (&now, TIME_UTC);
			printf ("\nRenversement de bits 4\nTemps écoulé :\t %lf ms\n",
			(now.tv_nsec - bf.tv_nsec) / 1e3);
  
  			timespec_get (&bf, TIME_UTC);
  			fftf_rdx4 (sig, TF, N, twiddles,cmp_add,cmp_mul);
			timespec_get (&now, TIME_UTC);
			printf ("\n");
			printf ("\n\tRADIX 4\nTemps écoulé :\t %*lf ms\n", 4,
				  (now.tv_nsec - bf.tv_nsec) / 1e3);
			printf("\n%d additions complexes\n%d multiplications complexes\n",*cmp_add,*cmp_mul );
			}
			*cmp_mul = 0;
			*cmp_add= 0 ;

			state = !state;
  			free (TF);
  			free (twiddles);
  			free (sig);
		}
	}else{
		
		if(strspn(argv[1], "0123456789") == strlen(argv[1])) {
			tmp = strtol(argv[1],NULL,10);
			N = 1<<(tmp);
			twiddles = get_twiddles (N);
                        sig = calloc (N, sizeof (short) * (N));
                        TF =   (float *) calloc (N<<1, sizeof (float) * (N<<1));
                        if (!sig || !twiddles || !TF|| !cmp_add || !cmp_mul)
                        {
                                printf ("Pas assez de mémoire\n");
                              return 0;
                        }

                        printf ("\n\tTF sur %d termes ",N);

                        printf ("\nBesoins mémoires :");
                        printf ("\nTF :\t%d o \n", N * 8);
                        printf ("TWI : \t%d o \n", N * 8); 
                        printf ("SIG : \t%d o \n", N * 2); 
                        for (i = 0; i < N; i++)
                            *(sig + i) = i; 

                        timespec_get (&bf, TIME_UTC);
                        rvs_16_rdx2 (sig,N);
                        timespec_get (&now, TIME_UTC);
                        printf ("\nRenversement de bits 2\nTemps écoulé :\t %lf ms\n",
                        (now.tv_nsec - bf.tv_nsec) / 1e3);

                        timespec_get (&bf, TIME_UTC);
                        fftf_rdx2 (sig, TF, N, twiddles,cmp_add,cmp_mul);
                        timespec_get (&now, TIME_UTC);
                        printf ("\n\tRADIX 2\nTemps écoulé :\t %lf ms",
                        (now.tv_nsec - bf.tv_nsec) / 1e3);
                        printf("\n%d additions complexes\n%d multiplications complexes\n",*cmp_add,*cmp_mul );

		if(!(atoi(argv[1])&0x1)){
			*cmp_mul = 0;
			*cmp_add= 0 ;
  			timespec_get (&bf, TIME_UTC);
			rvs_16_rdx4 (sig, N);
			timespec_get (&now, TIME_UTC);
			printf ("\nRenversement de bits 4\nTemps écoulé :\t %lf ms\n",
			(now.tv_nsec - bf.tv_nsec) / 1e3);
  
  			timespec_get (&bf, TIME_UTC);
  			fftf_rdx4 (sig, TF, N, twiddles,cmp_add,cmp_mul);
			timespec_get (&now, TIME_UTC);
			printf ("\n");
			printf ("\n\tRADIX 4\nTemps écoulé :\t %*lf ms\n", 4,
				  (now.tv_nsec - bf.tv_nsec) / 1e3);
			printf("\n%d additions complexes\n%d multiplications complexes\n",*cmp_add,*cmp_mul );
		
		}


  			free (TF);
  			free (twiddles);
  			free (sig);
	    }
	
	}
	free(cmp_add);
	free(cmp_mul);
  return 0;
}
