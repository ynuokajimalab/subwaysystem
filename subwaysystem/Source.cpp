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

int main(void)
{
	STEREO_PCM pcm0, pcm1;
	int n, k, J, N, offset, frame, number_of_frame;
	double fe1, fe2 ,delta, threshold, *b, *w, *x_real, *x_imag, *A, *T, *y_real, *y_imag;

	stereo_wave_read(&pcm0, "thelatest.WAV"); /* WAVE�t�@�C������X�e���I�̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	threshold = 0.0; /* �������l */
	fe1 = 1000.0 / (pcm0.fs * 2);
	fe2 = 2000.0 / (pcm0.fs * 2);
	delta = 900.0 / (pcm0.fs * 2);

	J = (int)(3.1 / delta + 0.5) - 1;
	if (J % 2 == 1)
	{
		J++;
	}

	b = (double*)calloc((J + 1), sizeof(double));
	w = (double*)calloc((J + 1), sizeof(double));

	Hanning_window(w, (J + 1));

	FIR_BPF(fe1, fe2, J, b, w);

	N = 1024; /* DFT�̃T�C�Y */

	number_of_frame = (pcm0.length * 2 - N / 2) / (N / 2); /* �t���[���̐� */
	w = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	A = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	T = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */

	Hanning_window(w, N); /* �n�j���O�� */

	stereo_wave_read(&pcm0, "test2.WAV"); /* WAVE�t�@�C������X�e���I�̉��f�[�^����͂��� */

	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = N / 2 * frame;

		/* x(n)��FFT */
		for (n = 0; n < N / 2; n++)
		{
			x_real[2 * n] = pcm0.sL[offset / 2 + n] * w[n];
			x_real[2 * n + 1] = pcm0.sR[offset / 2 + n] * w[n];
			x_imag[2 * n] = 0.0;
			x_imag[2 * n + 1] = 0.0;
		}
		FFT(x_real, x_imag, N);

		/* �U���X�y�N�g���ƈʑ��X�y�N�g�� */
		for (k = 0; k < N; k++)
		{
			A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
			if (x_imag[k] != 0.0 && x_real[k] != 0.0)
			{
				T[k] = atan2(x_imag[k], x_real[k]);
			}
		}
		
		/* �X�y�N�g���T�u�g���N�V���� */
		for (k = 0; k < N; k++)
		{
			y_real[k] = A[k] * cos(T[k]);
			y_imag[k] = A[k] * sin(T[k]);
		}
		IFFT(y_real, y_imag, N);
		if (A[k] = 1.0)
			printf("time=%d", k);

		/* ���H���ʂ̘A�� */
		for (n = 0; n < N / 2; n++)
		{
			pcm1.sL[offset / 2 + n] += y_real[2 * n];
			pcm1.sR[offset / 2 + n] += y_real[2 * n + 1];
		}
	}
	if (A[k] = 1.0)
		printf("time=%d", k);

	stereo_wave_write(&pcm1, "try.WAV"); /* WAVE�t�@�C���ɃX�e���I�̉��f�[�^���o�͂��� */

	free(pcm0.sL); /* �������̉�� */
	free(pcm0.sR); /* �������̉�� */
	free(pcm1.sL); /* �������̉�� */
	free(pcm1.sR); /* �������̉�� */
	free(w); /* �������̉�� */
	free(x_real); /* �������̉�� */
	free(x_imag); /* �������̉�� */
	free(A); /* �������̉�� */
	free(T); /* �������̉�� */
	free(y_real); /* �������̉�� */
	free(y_imag); /* �������̉�� */

	return 0;
}

