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

#define N 500

int main(void)
{
	MONO_PCM pcm0, pcm1;
	int i, j, n, m, k, L, offset, frame, number_of_frame, number_fe1, number_fe2, count;
	double fe1, fe2, *w, *x_real, *x_imag, *y_real, *y_imag, *w_real, *w_imag, max, threshold, temp, sum, sum_total, sum_average, sum_all, alpha;
	double noise_time1, noise_time2, number_noise1, number_noise2;
	double data[N];
	double data_all[N];

	mono_wave_read(&pcm0, "ktsyk_snyk.wav"); /* WAVEファイルからモノラルの音データを入力する */

	pcm1.fs = pcm0.fs; /* 標本化周波数 */
	pcm1.bits = pcm0.bits; /* 量子化精度 */
	pcm1.length = pcm0.length; /* 音データの長さ */
	pcm1.s = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */

	fe1 = 1500.0 / pcm0.fs; /* エッジ周波数 */
	fe2 = 2000.0 / pcm0.fs; /* エッジ周波数 */
	L = 512; /* フレームの長さ */

	w = (double*)calloc(L, sizeof(double)); /* メモリの確保 */

	Hanning_window(w, L); /* ハニング窓 */


	number_of_frame = pcm0.length / L; /* フレームの数 */

	w_real = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	w_imag = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	x_real = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	x_imag = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	y_real = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	y_imag = (double*)calloc(L, sizeof(double)); /* メモリの確保 */

	count = 0; /* カウントされた回数 */
	number_fe1 = (int)floor(L * fe1 / 2);
	number_fe2 = (int)floor(L * fe2 / 2);

	noise_time1 = 0.0;
	noise_time2 = 7.0;
	number_noise1 = noise_time1 * pcm0.fs / (2 * L);
	number_noise2 = noise_time2 * pcm0.fs / (2 * L);
	sum = 0.0; /* 初期化 */
	sum_total = 0.0;
	alpha = 1.5;

	for (frame = number_noise1; frame < number_noise2 + 1; frame++)
	{
		offset = L * frame;

		/* x(n)のFFT */
		for (n = 0; n < L; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
		}
		for (n = 0; n < L; n++)
		{
			x_real[n] = pcm0.s[offset + n];
		}

		/* w(m)のFFT */
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

		/* 掛け合わせ */
		for (k = 0; k < L; k++)
		{
			y_real[k] = x_real[k] * w_real[k];
			y_imag[k] = x_imag[k] * w_imag[k];
		}
		FFT(y_real, y_imag, L);

		/* 上1/4の合計を利用した判定 */
		for (i = number_fe1; i < number_fe2 + 1; i++)
		{
			data[i] = y_real[i];
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

		/* 最大値を利用した判定 */
		/*max = y_real[number_fe1];
		for (i = number_fe1; i < number_fe2 + 1; i++)
		{
			if (y_real[i] > max)
			{
				max = y_real[i];
			}
		}*/
	}
	sum_total += sum;
	sum_average = sum_total / (number_noise2 - number_noise1 + 1);
	threshold = alpha * sum_average;

	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;

			/* x(n)のFFT */
			for (n = 0; n < L; n++)
			{
				x_real[n] = 0.0;
				x_imag[n] = 0.0;
			}
			for (n = 0; n < L; n++)
			{
				x_real[n] = pcm0.s[offset + n];
			}
			
			/* w(m)のFFT */
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

			/* 掛け合わせ */
			for (k = 0; k < L; k++)
			{
				y_real[k] = x_real[k] * w_real[k];
				y_imag[k] = x_imag[k] * w_imag[k];
			}
			FFT(y_real, y_imag, L);

			/* 上位半分の合計を利用した判定 */
			for (i = number_fe1; i < number_fe2 + 1; i++)
			{
				data_all[i] = y_real[i];
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

			sum_all = 0; /* 初期化 */
			for (i = number_fe1; i < (number_fe2 - number_fe1 + 1) / 4 + number_fe1 + 1; i++)
			{
				sum_all += data_all[i];
			}

			if (sum_all > threshold)
			{
				count++;
				printf("time:%d\n", offset * 2 /pcm0.fs);
			}

			/* 最大値を利用した判定 */
			//max = y_real[number_fe1];
			//for (i = number_fe1; i < number_fe2 + 1; i++)
			//{
			//	if (y_real[i] > max)
			//	{
			//		max = y_real[i];
			//	}
			//}

			//threshold = 16.0; /* しきい値 */

			//if (max > threshold)
			//{
			//	count++;
			//	printf("frame:%d	max = %lf\n", frame,max);
			//}

		//		/* フィルタリング結果の連結 */
		//		for (n = 0; n < L * 2; n++)
		//		{
		//			if (offset + n < pcm1.length)
		//			{
		//				pcm1.s[offset + n] += y_real[n];
		//			}
		//		}
		//	}
	}
	printf("回数は%d回です\n", count);

	//mono_wave_write(&pcm1, "short1_out.WAV"); /* WAVEファイルにモノラルの音データを出力する */

	free(pcm0.s); /* メモリの解放 */
	free(pcm1.s); /* メモリの解放 */
	free(w); /* メモリの解放 */
	free(w_real); /* メモリの解放 */
	free(w_imag); /* メモリの解放 */
	free(x_real); /* メモリの解放 */
	free(x_imag); /* メモリの解放 */
	free(y_real); /* メモリの解放 */
	free(y_imag); /* メモリの解放 */

	return 0;
}
