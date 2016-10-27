#pragma warning(disable:4996)

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "wave.h"
#include "window_function.h"
#include "sinc.h"
#include "fir_filter.h"
#include "fft.h"
#include "countSound.h"

//int main(void)
//{
//	STEREO_PCM pcm0, pcm1;
//	int n, m, k, J, L, N, offset, frame, number_of_frame,countsound;
//	double fe, delta, threshold, *b, *w, *b_real, *b_imag, *x_real, *x_imag, *y_real, *y_imag;
//
//	stereo_wave_read(&pcm0, "005_160615_0941V0.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */
//
//	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
//	pcm1.bits = pcm0.bits; /* �ʎq�����x */
//	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
//	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�ہA���� */
//	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� �A����*/
//
//	fe = 500.0 / pcm0.fs; /* �G�b�W���g�� */
//	delta = 600.0 / pcm0.fs; /* �J�ڑш敝 */
//
//	J = (int)(3.1 / delta + 0.5) - 1; /* �x����̐� */
//	if (J % 2 == 1)
//	{
//		J++; /* J+1����ɂȂ�悤�ɒ������� */
//	}
//
//	b = (double*) calloc((J + 1), sizeof(double)); /* �������̊m�� */
//	w = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */
//
//	Hanning_window(w, (J + 1)); /* �n�j���O�� */
//
//	FIR_LPF(fe, J, b, w); /* FIR�t�B���^�̐݌v */
//
//	L = 256; /* �t���[���̒��� */
//	N = 512; /* DFT�̃T�C�Y */
//
//	countsound = 0; /* ����臒l�𒴂�����*/
//	threshold = 100.0;/* 臒l*/
//
//	number_of_frame = pcm0.length / L; /* �t���[���̐� */
//
//	b_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
//	b_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
//	x_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
//	x_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
//	y_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
//	y_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
//
//	/*for (m = 0; m < J; m++)
//		{
//
//			for (n = 0; n < N; n++)
//			{
//				y_real[n] += b_real[m] * x_real[n - m];
//				y_imag[n] += b_imag[m] * x_imag[n - m];
//			}
//
//		}*/
//	 
//	  for (frame = 0; frame < number_of_frame; frame++)
//	  {
//	  	offset = L * frame;
//
//
//	  	/* x(n)��FFT */
//	  	  /*for (n = 0; n < N; n++)
//	  	{
//	  		x_real[n] = 0.0;
//	  		x_imag[n] = 0.0;
//		}
//		for (n = 0; n < L/2; n++)
//	  	{
//	  		x_real[2 * n] = pcm0.sL[offset/2 + n];
//	  		x_real[2 * n + 1] = pcm0.sR[offset/2 + n];
//	  	}*/
//	  	/*FFT(x_real, x_imag, N);*/
//
//	  	/* b(m)��FFT */
//	  	  for (m = 0; m < N; m++)
//	  	{
//	  		b_real[m] = 0.0;
//	  		b_imag[m] = 0.0;
//	  	}
//	  	for (m = 0; m <= J; m++)
//	  	{
//	  		b_real[m] = b[m];
//	  	}
//	  	/*FFT(b_real, b_imag, N);*/
//
//	  	countsound += judgeSounnd(y_real, N, threshold);
//
//	  	/* �t�B���^�����O */
//		for (n = 0; n < pcm0.length; n++)
//		{
//			for (m = 0; m <= J; m++)
//			{
//				if (n - m >= 0)
//				{
//					pcm1.sL[n] += b[m] * pcm0.sL[n - m];
//					pcm1.sR[n] += b[m] * pcm0.sR[n - m];
//				}
//			}
//		}
//	  
//
//	  	
//	  	
//	  	/*IFFT(y_real, y_imag, N);*/
//
//	  	/* �t�B���^�����O���ʂ̘A�� */
//	  	for (n = 0; n < L; n++)
//	  	{
//	  		if (offset/2 + n < pcm1.length)
//	  		{
//	  			pcm1.sL[offset/2 + n] += y_real[2 * n];
//	  			pcm1.sR[offset/2 + n] += y_real[2 * n + 1];
//	  		}
//	  	}
//	  
//	  }
//	stereo_wave_write(&pcm1, "ex6_4extra.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
//
//	free(pcm0.sL); /* �������̉�� */
//	free(pcm0.sR); /* �������̉�� */
//	free(pcm1.sL); /* �������̉�� */
//	free(pcm1.sR); /* �������̉�� */
//	free(b); /* �������̉�� */
//	free(w); /* �������̉�� */
//	free(b_real); /* �������̉�� */
//	free(b_imag); /* �������̉�� */
//	free(x_real); /* �������̉�� */
//	free(x_imag); /* �������̉�� */
//	free(y_real); /* �������̉�� */
//	free(y_imag); /* �������̉�� */
//
//	printf("countsound=%d\n",countsound);
//	printf("finish");
//
//	return 0;
//}
int main(void)
{
	STEREO_PCM pcm0, pcm1;
	int n, m, J, L, offset, frame, number_of_frame;
	double fe, delta, *b, *w, *x, *y;

	stereo_wave_read(&pcm0, "005_160615_0941V0.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	fe = 450.0 / (pcm0.fs * 2); /* �G�b�W���g�� */
	delta = 450.0 / (pcm0.fs * 2); /* �J�ڑш敝 */

	J = (int)(3.1 / delta + 0.5) - 1; /* �x����̐� */
	if (J % 2 == 1)
	{
		J++; /* J+1����ɂȂ�悤�ɒ������� */
	}

	b = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */
	w = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */

	Hanning_window(w, (J + 1)); /* �n�j���O�� */

	FIR_LPF(fe, J, b, w); /* FIR�t�B���^�̐݌v */

	L = 256; /* �t���[���̒��� */

	number_of_frame = (pcm0.length * 2) / L; /* �t���[���̐� */

	x = (double*)calloc((L + J), sizeof(double)); /* �������̊m�� */
	y = (double*)calloc(L, sizeof(double)); /* �������̊m�� */

	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;

		/* ���O�̃t���[���̌㔼��J�T���v������������ */
		for (n = 0; n <(L + J) / 2; n++)
		{
			if ((offset / 2) - J / 2 + n < 0)
			{
				x[n] = 0.0;/* �ŏ��̃t���[���̒��O��J�T���v���͎��ۂɂ͑��݂��Ȃ� */
			}
			else
			{
				x[2 * n] = pcm0.sL[(offset / 2) - J / 2 + n];
				x[2 * n + 1] = pcm0.sR[(offset / 2) - J / 2 + n];
			}
		}

		for (n = 0; n < L; n++)
		{
			y[n] = 0.0;
		}

		/* �t�B���^�����O */
		for (n = 0; n < L; n++)
		{
			for (m = 0; m <= J; m++)
			{
				y[n] += b[m] * x[J + n - m];
			}
		}

		/* �t�B���^�����O���ʂ̘A�� */
		for (n = 0; n < L / 2; n++)
		{
			pcm1.sL[(offset / 2) + n] = y[2 * n];
			pcm1.sR[(offset / 2) + n] = y[(2 * n) + 1];
		}
	}

	stereo_wave_write(&pcm1, "ex6_4test.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */

	free(pcm0.sL); /* �������̉�� */
	free(pcm0.sR); /* �������̉�� */
	free(pcm1.sL); /* �������̉�� */
	free(pcm1.sR); /* �������̉�� */
	free(b); /* �������̉�� */
	free(w); /* �������̉�� */
	free(x); /* �������̉�� */
	free(y); /* �������̉�� */

	return 0;
}

