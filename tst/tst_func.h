/******************************************************************************
                        PROJET  S6
        EL-KHARROUBI    GRAVES  LEFEBVRE        NEZET
*******************************************************************************/

#ifndef __TST_FUNC__
#define __TST_FUNC__

/*!
 * @function    br_tst
 * @abstract    fonction permettant le test d'une fonction de renversement de bit
 * @param       br      la fonction de renversement
 *
 */
float
br_tst (void (*br) (short *, int), int N, short *sig, struct timespec *bf,
	struct timespec *now);
/*!
 * @function    fftf_tst 
 * @abstract    fonction permettant le test d'une fft
 * @param       fft     pointeur vers la fft
 */

float
fftf_tst (void (*fft)
	  (short *, float *, int, float *, unsigned int *, unsigned int *),
	  int N, short *sig, float *TF, float *twiddles, struct timespec *bf,
	  struct timespec *now, unsigned int *add, unsigned int *mul);
float
fftf_cmp_tst (void (*fft)
	      (short *, complex *, int, complex *, unsigned int *,
	       unsigned int *), int N, short *sig, complex * TF,
	      complex * twiddles, struct timespec *bf, struct timespec *now,
	      unsigned int *add, unsigned int *mul);
float
fftw_tst (int N, short *sig, float *TF, struct timespec *bf,
	  struct timespec *now);
int no_arg (struct timespec bf, struct timespec now, FILE * res);
int tst_one_tf (char *number, struct timespec bf, struct timespec now,
		FILE * res);

int init (char *argv2);
int ml_compare (char *sig_name, char *ml_name,FILE* res,struct timespec bf, struct timespec now);
#endif
