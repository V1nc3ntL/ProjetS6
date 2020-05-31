/******************************************************************************
			PROJET 	S6
	EL-KHARROUBI 	GRAVES 	LEFEBVRE 	NEZET
*******************************************************************************/
#ifndef __FFT_TST_H__
#define __FFT_TST_H__
void fftf_rdx2_tst (short *sig, float *TF, int N, float *twiddles,unsigned int * cmp_add, unsigned int * cmp_mul);
void fftf_rdx4_tst (short *sig, float *TF, int N, float *twiddles,unsigned int * cmp_add, unsigned int * cmp_mul);
void fftf_rdx2_tst_re_im (short *sig, float *TF, int N, float *twiddles,unsigned int * cmp_add, unsigned     int * cmp_mul);
float * get_twiddles (int b);
#endif
