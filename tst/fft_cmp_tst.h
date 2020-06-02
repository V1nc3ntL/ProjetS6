/******************************************************************************
			PROJET 	S6
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
*******************************************************************************/
#ifndef __FFT_CMP_TST_H__
#define __FFT_CMP_TST_H__
void fftf_rdx2_cmp_tst (short *sig, complex * TF, int N, complex * twiddles,
			unsigned int *cmp_add, unsigned int *cmp_mul);
void fftf_rdx4_cmp_tst (short *sig, complex * TF, int N, complex * twiddles,
			unsigned int *cmp_add, unsigned int *cmp_mul);
complex *get_twiddles_cmp (int b);
#endif
