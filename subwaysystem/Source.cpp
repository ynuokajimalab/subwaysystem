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

	mono_wave_read(&pcm0, "ktsy_syk_33.wav"); /* WAVEファイルからモノラルの音データを入力する */

	pcm1.fs = pcm0.fs; /* 標本化周波数 */
	pcm1.bits = pcm0.bits; /* 量子化精度 */
	pcm1.length = pcm0.length; /* 音データの長さ */
	pcm1.s = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */

	fe1 = 1500.0 / pcm0.fs; /* エッジ周波数 */
	fe2 = 2000.0 / pcm0.fs; /* エッジ周波数 */
	delta = 500.0 / pcm0.fs; /* 遷移帯域幅 */

	J = (int)(3.1 / delta + 0.5) - 1; /* 遅延器の数 */
	if (J % 2 == 1)
	{
		J++; /* J+1が奇数になるように調整する */
	}

	b = (double*)calloc((J + 1), sizeof(double)); /* メモリの確保 */
	w = (double*)calloc((J + 1), sizeof(double)); /* メモリの確保 */

	Hanning_window(w, (J + 1)); /* ハニング窓 */

	//FIR_BPF(fe1,fe2, J, b, w); /* FIRフィルタの設計 */

	L = 256; /* フレームの長さ */
	N = 512; /* DFTのサイズ */

	number_of_frame = pcm0.length / L; /* フレームの数 */

	w_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	w_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */

	count = 0; /* カウントされた回数 */
	number_fe1 = floor(N * fe1 / (pcm0.fs / 2));
	number_fe2 = floor(N * fe2 / (pcm0.fs / 2));

	printf("fe1 = %lf,fe2 = %lf\n",number_fe1,number_fe2);

	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;

			/* x(n)のFFT */
			for (n = 0; n < N; n++)
			{
				x_real[n] = 0.0;
				x_imag[n] = 0.0;
			}
			for (n = 0; n < L; n++)
			{
				x_real[n] = pcm0.s[offset + n];
			}
			
			/* w(m)のFFT */
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

			/* 掛け合わせ */
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

			threshold = 5.0; /* しきい値 */

			if (max > threshold)
			{
				count++;
			}

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

	free(pcm0.s); /* メモリの解放 */
	free(pcm1.s); /* メモリの解放 */
	free(b); /* メモリの解放 */
	free(w); /* メモリの解放 */
	free(w_real); /* メモリの解放 */
	free(w_imag); /* メモリの解放 */
	free(x_real); /* メモリの解放 */
	free(x_imag); /* メモリの解放 */
	free(y_real); /* メモリの解放 */
	free(y_imag); /* メモリの解放 */

	return 0;
}
