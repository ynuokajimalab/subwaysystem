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
	MONO_PCM pcm0, pcm1;
	int n, m, k, J, L, N, offset, frame, number_of_frame,countsound;
	double fe, delta, threshold, *b, *w, *b_real, *b_imag, *x_real, *x_imag, *y_real, *y_imag;

	mono_wave_read(&pcm0, "005_160615_0941V0.wav"); /* WAVEファイルからモノラルの音データを入力する */

	pcm1.fs = pcm0.fs; /* 標本化周波数 */
	pcm1.bits = pcm0.bits; /* 量子化精度 */
	pcm1.length = pcm0.length; /* 音データの長さ */
	pcm1.s = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */

	fe = 600.0 / pcm0.fs; /* エッジ周波数 */
	delta = 500.0 / pcm0.fs; /* 遷移帯域幅 */

	J = (int)(3.1 / delta + 0.5) - 1; /* 遅延器の数 */
	if (J % 2 == 1)
	{
		J++; /* J+1が奇数になるように調整する */
	}

	b = (double*) calloc((J + 1), sizeof(double)); /* メモリの確保 */
	w = (double*)calloc((J + 1), sizeof(double)); /* メモリの確保 */

	Hanning_window(w, (J + 1)); /* ハニング窓 */

	FIR_LPF(fe, J, b, w); /* FIRフィルタの設計 */

	L = 256; /* フレームの長さ */
	N = 512; /* DFTのサイズ */

	countsound = 0; /* 音が閾値を超えた回数*/
	threshold = 100.0;/* 閾値*/

	number_of_frame = pcm0.length / L; /* フレームの数 */

	b_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	b_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */

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
		FFT(x_real, x_imag, N);

		/* b(m)のFFT */
		for (m = 0; m < N; m++)
		{
			b_real[m] = 0.0;
			b_imag[m] = 0.0;
		}
		for (m = 0; m <= J; m++)
		{
			b_real[m] = b[m];
		}
		FFT(b_real, b_imag, N);

		/* フィルタリング */
		for (k = 0; k < N; k++)
		{
			y_real[k] = x_real[k] * b_real[k] - x_imag[k] * b_imag[k];
			y_imag[k] = x_imag[k] * b_real[k] + x_real[k] * b_imag[k];
		}

		countsound += judgeSounnd(y_real,N,threshold);
		
		IFFT(y_real, y_imag, N);

		/* フィルタリング結果の連結 */
		for (n = 0; n < L * 2; n++)
		{
			if (offset + n < pcm1.length)
			{
				pcm1.s[offset + n] += y_real[n];
			}
		}
	}

	mono_wave_write(&pcm1, "ex6_4.wav"); /* WAVEファイルにモノラルの音データを出力する */

	free(pcm0.s); /* メモリの解放 */
	free(pcm1.s); /* メモリの解放 */
	free(b); /* メモリの解放 */
	free(w); /* メモリの解放 */
	free(b_real); /* メモリの解放 */
	free(b_imag); /* メモリの解放 */
	free(x_real); /* メモリの解放 */
	free(x_imag); /* メモリの解放 */
	free(y_real); /* メモリの解放 */
	free(y_imag); /* メモリの解放 */

	printf("countsound =%d\n",countsound);
	printf("finish");

	return 0;
}
