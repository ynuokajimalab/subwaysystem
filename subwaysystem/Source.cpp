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

#define N 500

int main(void)
{
	MONO_PCM pcm0, pcm1;
	int i, j, n, m, k, L, offset, frame, number_of_frame, number_fe1, number_fe2, count;
	double fe1, fe2, *w, *x_real, *x_imag, *y_real, *y_imag, *w_real, *w_imag, max, max_total, max_average, threshold, temp, sum, sum_total, sum_average, sum_all, alpha, *A, *A_all, max_top;
	double noise_time1, noise_time2, number_noise1, number_noise2;
	double data[N];
	double data_all[N];

	mono_wave_read(&pcm0, "30_29.wav"); /* WAVE�t�@�C�����烂�m�����̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.s = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	/* ���g���͈̔͂��w�� */
	fe1 = 500.0 / pcm0.fs;
	fe2 = 900.0 / pcm0.fs;
	L = 256; /* �t���[���̒��� */

	w = (double*)calloc(L, sizeof(double)); /* �������̊m�� */

	Hanning_window(w, L); /* �n�j���O�� */


	number_of_frame = pcm0.length / L; /* �t���[���̐� */
	printf("number_of_frame:%d\n", number_of_frame);

	w_real = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	w_imag = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	x_real = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	x_imag = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	y_real = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	y_imag = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	w_real = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	w_imag = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	A = (double*)calloc(L, sizeof(double)); /* �������̊m�� */
	A_all = (double*)calloc(L, sizeof(double)); /* �������̊m�� */

	count = 0; /* �J�E���g���ꂽ�� */
	number_fe1 = (int)floor(L * fe1);
	number_fe2 = (int)floor(L * fe2);

	printf("number_fe1:%d     number_fe2:%d\n", number_fe1,number_fe2);

	noise_time1 = 3.9;
	noise_time2 = 4.1;
	number_noise1 = noise_time1 * pcm0.fs /  L;
	number_noise2 = noise_time2 * pcm0.fs /  L;
	sum = 0.0; /* ������ */
	sum_total = 0.0;
	alpha = 1.1;
	max_total = 0.0;

	for (frame = number_noise1; frame < number_noise2 + 1; frame++)
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

		for (k = 0; k < L; k++)
		{
			A[k] = sqrt(y_real[k] * y_real[k] + y_imag[k] * y_imag[k]);
		}

		/* ��1/4�̍��v�𗘗p�������� */
		for (i = number_fe1; i < number_fe2 + 1; i++)
		{
			data[i] = A[i];
			for (j = i + 1; j < number_fe2 + 1; j++)
			{
				if (data[i] < data[j])
				{
					temp = data[i];
					data[i] = data[j];
					data[j] = temp;
			    }
			}
		}
		for (i = number_fe1; i < (number_fe2 - number_fe1 + 1) / 4 + number_fe1 + 1; i++)
		{
			sum += data[i];
		}
		IFFT(y_real, y_imag, L);

		/* �ő�l�𗘗p�������� */
		/*for (i = number_fe1; i < number_fe2 + 1; i++)
		{
			data[i] = A[i];
        }
		max = y_real[number_fe1];
		for (i = number_fe1; i < number_fe2 + 1; i++)
		{
			if (data[i] > max)
			{
				max = data[i];
			}
		}
		max_top = max;
		if (max > max_top)
		{
			max_top = max;
		}*/
	}
	sum_total += sum;
	sum_average = sum_total / (number_noise2 - number_noise1 + 1);
	threshold = alpha * sum_average;

	/*max_total += max;
	max_average = max_total / (number_noise2 - number_noise1 + 1);
	threshold = alpha * max_top;*/

	for (frame = 0; frame < number_of_frame/2; frame++)
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

			for (k = 0; k < L; k++)
			{
				A_all[k] = sqrt(y_real[k] * y_real[k] + y_imag[k] * y_imag[k]);
			}

			/* ���1/4�̍��v�𗘗p�������� */
			for (i = number_fe1; i < number_fe2 + 1; i++)
			{
				data_all[i] = A_all[i];
				for (j = i + 1; j < number_fe2 + 1; j++)
				{
					if (data_all[i] < data_all[j])
					{
						temp = data_all[i];
						data_all[i] = data_all[j];
						data_all[j] = temp;
					}
				}
			}

			sum_all = 0; /* ������ */
			for (i = number_fe1; i < (number_fe2 - number_fe1 + 1) / 4 + number_fe1 + 1; i++)
			{
				sum_all += data_all[i];
			}

			if (sum_all > threshold)
			{
				count++;
				//?
				printf("time:%f\n", (double)offset*2 /pcm0.fs);
			}

			/* �ő�l�𗘗p�������� */
			/*for (i = number_fe1; i < number_fe2 + 1; i++)
			{
				data_all[i] = A_all[i];
			}

			max = data_all[number_fe1];
			for (i = number_fe1; i < number_fe2 + 1; i++)
			{
				if (data_all[i] > max)
				{
					max = data_all[i];
				}
			}

			if (max > threshold)
			{
				count++;
				printf("time:%f\n", (double)offset * 2 / pcm0.fs);
			}*/

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

	mono_wave_write(&pcm1, "test1.WAV"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */

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
