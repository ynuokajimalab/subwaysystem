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

int main(void)
{
	MONO_PCM pcm0, pcm1;
	int i, n, m, k, J, L, N, offset, frame, number_of_frame, number_fe1, number_fe2, count;
	double fe1, fe2, delta, *b, *w, *b_real, *b_imag, *x_real, *x_imag, *y_real, *y_imag, *w_real, *w_imag, max, threshold;

	mono_wave_read(&pcm0, "ktsy_syk_33.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.s = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	fe1 = 1500.0 / pcm0.fs; /* �G�b�W���g�� */
	fe2 = 2000.0 / pcm0.fs; /* �G�b�W���g�� */
	delta = 500.0 / pcm0.fs; /* �J�ڑш敝 */

	J = (int)(3.1 / delta + 0.5) - 1; /* �x����̐� */
	if (J % 2 == 1)
	{
		J++; /* J+1����ɂȂ�悤�ɒ������� */
	}

	b = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */
	w = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */

	Hanning_window(w, (J + 1)); /* �n�j���O�� */

	//FIR_BPF(fe1,fe2, J, b, w); /* FIR�t�B���^�̐݌v */

	L = 256; /* �t���[���̒��� */
	N = 512; /* DFT�̃T�C�Y */

	number_of_frame = pcm0.length / L; /* �t���[���̐� */

	w_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	w_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */

	count = 0; /* �J�E���g���ꂽ�� */
	number_fe1 = floor(N * fe1 / (pcm0.fs / 2));
	number_fe2 = floor(N * fe2 / (pcm0.fs / 2));

	printf("fe1 = %lf,fe2 = %lf\n",number_fe1,number_fe2);

	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;

			/* x(n)��FFT */
			for (n = 0; n < N; n++)
			{
				x_real[n] = 0.0;
				x_imag[n] = 0.0;
			}
			for (n = 0; n < L; n++)
			{
				x_real[n] = pcm0.s[offset + n];
			}
			
			/* w(m)��FFT */
			for (m = 0; m < N; m++)
			{
				w_real[m] = 0.0;
				w_imag[m] = 0.0;
			}
			for (m = 0; m <= J; m++)
			{
				w_real[m] = w[m];
			}
			/*FFT(b_real, b_imag, N);*/

			/* �|�����킹 */
			for (k = 0; k < N; k++)
			{
				y_real[i] = x_real[k] * w_real[k];
				y_imag[i] = x_imag[k] * w_imag[k];
			}
			FFT(y_real, y_imag, N);

			max = y_real[number_fe1];
			for (i = number_fe1; i < number_fe2 + 1; i++)
			{
				if (y_real[i] > max)
				{
					max = y_real[i];
				}
			}

			threshold = 5.0; /* �������l */

			if (max > threshold)
			{
				count++;
			}

		//		/* �t�B���^�����O���ʂ̘A�� */
		//		for (n = 0; n < L * 2; n++)
		//		{
		//			if (offset + n < pcm1.length)
		//			{
		//				pcm1.s[offset + n] += y_real[n];
		//			}
		//		}
		//	}
	}
	printf("�񐔂�%d��ł�\n", count);

	free(pcm0.s); /* �������̉�� */
	free(pcm1.s); /* �������̉�� */
	free(b); /* �������̉�� */
	free(w); /* �������̉�� */
	free(w_real); /* �������̉�� */
	free(w_imag); /* �������̉�� */
	free(x_real); /* �������̉�� */
	free(x_imag); /* �������̉�� */
	free(y_real); /* �������̉�� */
	free(y_imag); /* �������̉�� */

	return 0;
}
