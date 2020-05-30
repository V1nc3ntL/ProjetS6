#ifndef __FFT_H__
#define __FFT_H__
void fftf_rdx2 (short *sig, float *TF, int N, float *twiddles);
void fftf_rdx4 (short *sig, float *TF, int N, float *twiddles);
float * get_twiddles (int b);
#endif
