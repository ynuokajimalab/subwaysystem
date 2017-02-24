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
	int i, j, n, m, k, L, offset, frame, number_of_frame, number_fe1, number_fe2, count;
	double fe1, fe2, *w, *x_real, *x_imag, *y_real, *y_imag, *w_real, *w_imag, max, max_top, threshold, temp, sum, max_top;
	double noise_time1, noise_time2, frame_noise1, frame_noise2;

	mono_wave_read(&pcm0, "minami_nakam.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.s = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	printf("fs = %d\n",pcm0.fs);

	fe1 = 1500.0 / pcm0.fs; /* �G�b�W���g�� */
	fe2 = 2000.0 / pcm0.fs; /* �G�b�W���g�� */
	L = 512; /* �t���[���̒��� */

	w = (double*)calloc(L, sizeof(double)); /* �������̊m�� */

	Hanning_window(w, L); /* �n�j���O�� */


	number_of_frame = pcm0.length / L; /* �t���[���̐� */

	w_real = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	w_imag = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	x_real = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	x_imag = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	y_real = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	y_imag = (double*)calloc(L, sizeof(double)); /* �������̊m�� */

	count = 0; /* �J�E���g���ꂽ�� */
	number_fe1 = (int)floor(L * fe1 / 2);
	number_fe2 = (int)floor(L * fe2 / 2);

	noise_time1 = 0.0;
	noise_time2 = 8.0;
	frame_noise1 = floor(noise_time1 * pcm0.fs / L);
	frame_noise2 = floor(noise_time2 * pcm0.fs / L);

	for (frame = frame_noise1; frame < frame_noise2 + 1; frame++)
	{
		offset = L * frame;

		/* x(n)��FFT */
		for (n = 0; n < L; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
		}
		for (n = 0; n < L; n++)
		{
			x_real[n] = pcm0.s[offset + n];
		}

		/* w(m)��FFT */
		for (m = 0; m < L; m++)
		{
			w_real[m] = 0.0;
			w_imag[m] = 0.0;
		}
		for (m = 0; m < L; m++)
		{
			w_real[m] = w[m];
		}
		/*FFT(b_real, b_imag, N);*/

		/* �|�����킹 */
		for (k = 0; k < L; k++)
		{
			y_real[k] = x_real[k] * w_real[k];
			y_imag[k] = x_imag[k] * w_imag[k];
		}
		FFT(y_real, y_imag, L);

		/* ��ʔ����̍��v�𗘗p�������� */

		max = y_real[number_fe1];

		for (j = number_fe1; j < number_fe2 + 1; j++)
		{
			if (y_real[j] >max)
			{
				max = y_real[j];
			}
		}
		max_top = max;

	}


	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;
		printf("time:%f\n", ((double)offset*2 / (double)pcm0.fs));

			/* x(n)��FFT */
			for (n = 0; n < L; n++)
			{
				x_real[n] = 0.0;
				x_imag[n] = 0.0;
			}
			for (n = 0; n < L; n++)
			{
				x_real[n] = pcm0.s[offset + n];
			}
			
			/* w(m)��FFT */
			for (m = 0; m < L; m++)
			{
				w_real[m] = 0.0;
				w_imag[m] = 0.0;
			}
			for (m = 0; m < L; m++)
			{
				w_real[m] = w[m];
			}
			/*FFT(b_real, b_imag, N);*/

			/* �|�����킹 */
			for (k = 0; k < L; k++)
			{
				y_real[k] = x_real[k] * w_real[k];
				y_imag[k] = x_imag[k] * w_imag[k];
			}
			FFT(y_real, y_imag, L);

			/* ��ʔ����̍��v�𗘗p�������� */
			for (j = number_fe1; j < number_fe2 + 1; j++)
			{
				if (y_real[j] < y_real[j + 1])
				{
					temp = y_real[j];
					y_real[j] = y_real[j + 1];
					y_real[j + 1] = temp;
				}
			}

			sum = 0; /* ������ */

			for (i = number_fe1; i < (number_fe2 - number_fe1 + 1) / 4 + number_fe1 + 1; i++)
			{
				sum += y_real[i];
			}
			if (sum > threshold)
			{
				count++;
				printf("time:%f\n", ((double)offset/(double)pcm0.fs));
			}

			///* �ő�l�𗘗p�������� */
			//max = y_real[number_fe1];
			//for (i = number_fe1; i < number_fe2 + 1; i++)
			//{
			//	if (y_real[i] > max)
			//	{
			//		max = y_real[i];
			//	}
			//}

			//threshold = 57.0; /* �������l */

			//if (max > threshold)
			//{
			//	count++;
			//	printf("frame:%d	max = %lf\n", frame,max);
			//}

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

	mono_wave_write(&pcm1, "test.WAV"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */

	free(pcm0.s); /* �������̉�� */
	free(pcm1.s); /* �������̉�� */
	free(w); /* �������̉�� */
	free(w_real); /* �������̉�� */
	free(w_imag); /* �������̉�� */
	free(x_real); /* �������̉�� */
	free(x_imag); /* �������̉�� */
	free(y_real); /* �������̉�� */
	free(y_imag); /* �������̉�� */

	return 0;
}
