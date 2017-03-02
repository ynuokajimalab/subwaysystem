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
#include "frame-time.h"

#define N 2000

int main(void)
{
	MONO_PCM pcm0, pcm1;
	int i, j, n, m, k, L, offset, frame, number_of_frame, number_fe1, number_fe2, count, n_noise1, n_noise2;
	double fe1, fe2, *w, *x_real, *x_imag, *y_real, *y_imag, *w_real, *w_imag, *A, *A_all, *xw_real, *xw_imag, *xw;
    double max, max_total, max_average, threshold, temp, sum, sum_total, sum_average, sum_all, alpha, max_top;
	double noise_time1, noise_time2, number_noise1, number_noise2;
	double data[N];

	mono_wave_read(&pcm0,"ftk_yg_16bit.wav"); /* WAVEファイルからモノラルの音データを入力する */
			
	pcm1.fs = pcm0.fs; /* 標本化周波数 */
	pcm1.bits = pcm0.bits; /* 量子化精度 */
	pcm1.length = pcm0.length; /* 音データの長さ */
	pcm1.s = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */

	for (int i = 0; i < pcm0.length; i++)
	{
		pcm1.s[i] = pcm0.s[i];
	}
	/* 周波数の範囲を指定 */
	fe1 = 0.0 / pcm0.fs;
	fe2 = 3000.0 / pcm0.fs;
	L = 1024; /* フレームの長さ */
	alpha = 40.0;

	w = (double*)calloc(L, sizeof(double)); /* メモリの確保 */

	Hanning_window(w, L); /* ハニング窓 */

	number_of_frame = pcm0.length / L; /* フレームの数 */
	printf("number_of_frame:%d\n", number_of_frame);

	w_real = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	w_imag = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	x_real = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	x_imag = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	y_real = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	y_imag = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	w_real = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	w_imag = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	xw_real = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	xw_imag = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	xw = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	A = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	A_all = (double*)calloc(L, sizeof(double)); /* メモリの確保 */

	count = 0; /* カウントされた回数 */
	number_fe1 = (int)floor(L * fe1);
	number_fe2 = (int)floor(L * fe2);

	noise_time1 = 0.0;
	noise_time2 = 0.2;
	number_noise1 = noise_time1 * pcm0.fs /  L;
	number_noise2 = noise_time2 * pcm0.fs /  L;
	n_noise1 = noise_time1 * pcm0.fs;
	n_noise2 = noise_time2 * pcm0.fs;
	
	sum = 0.0; /* 初期化 */
	sum_total = 0.0;
	max_total = 0.0;
	
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

		for (k = 0; k < L; k++)
		{
			A[k] = sqrt(y_real[k] * y_real[k] + y_imag[k] * y_imag[k]);
		}

		/* 上1/4の合計を利用した判定 */
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
		for (i = number_fe1; i < (number_fe2 - number_fe1) / 4 + number_fe1 + 1; i++)
		{
			sum += data[i];
		}
		IFFT(y_real, y_imag, L);

	}
	sum_total += sum;
	sum_average = sum_total / (number_noise2 - number_noise1 + 1);
	threshold = alpha * sum_average;
	printf("threshold:%f\n", threshold);
	
	CountTimes_mono(L, 0., 0.0, 0.0, number_fe1, number_fe2, threshold, 0.3, 2.5, &pcm0);

	mono_wave_write(&pcm1, "test1.WAV"); /* WAVEファイルにモノラルの音データを出力する */

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
