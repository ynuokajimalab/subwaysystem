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
	int m, n, J, k, N ,L, offset, offset_ex, frame, number_of_frame, number_of_frame_ex;
	double fe1, fe2, delta, threshold, *b, *w, *x, *x_real, *x_imag, *A, *T, *y, *y_real, *y_imag;

	stereo_wave_read(&pcm0, "thelatest.WAV"); /* WAVEファイルからステレオの音データを入力する */

	pcm1.fs = pcm0.fs; /* 標本化周波数 */
	pcm1.bits = pcm0.bits; /* 量子化精度 */
	pcm1.length = pcm0.length; /* 音データの長さ */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */

	fe1 = 1000.0 / (pcm0.fs * 2); /* エッジ周波数 */
    fe2 = 2000.0 / (pcm0.fs * 2); /* エッジ周波数 */
	delta = 450.0 / (pcm0.fs * 2); /* 遷移帯域幅 */

	J = (int)(3.1 / delta + 0.5) - 1; /* 遅延器の数 */
	if (J % 2 == 1)
	{
	        J++; /* J+1が奇数になるように調整する */
	}

	threshold = 0.5; /* しきい値 */

	L = 256; /* フレームの長さ */
	N = 1024; /* DFTのサイズ */

	number_of_frame = (pcm0.length * 2 - N / 2) / (N / 2); /* フレームの数 */
	number_of_frame_ex = (pcm0.length * 2) / L; /* フレームの数 */

	b = (double*)calloc((J + 1), sizeof(double)); /* メモリの確保 */
	w = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x = (double*)calloc((L + J), sizeof(double)); /* メモリの確保 */
	x_real =(double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	A = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	T = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y = (double*)calloc(L, sizeof(double)); /* メモリの確保 */
	y_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */

	Hanning_window(w, (J+1)); /* ハニング窓 */

	FIR_BPF(fe1, fe2, J, b, w); /* FIRフィルタの設計 */

	for (frame = 0; frame < number_of_frame_ex; frame++)
	{
		offset_ex = L * frame;

		/* 直前のフレームの後半のJサンプルをつけ加える */
		for (n = 0; n <(L + J) / 2; n++)
		{
			if ((offset_ex / 2) - J / 2 + n < 0)
			{
				x[n] = 0.0;/* 最初のフレームの直前のJサンプルは実際には存在しない */
			}
			else
			{
				x[2 * n] = pcm0.sL[(offset_ex / 2) - J / 2 + n];
				x[2 * n + 1] = pcm0.sR[(offset_ex / 2) - J / 2 + n];
			}
		}

		for (n = 0; n < L; n++)
		{
			y[n] = 0.0;
		}

		/* フィルタリング */
		for (n = 0; n < L; n++)
		{
			for (m = 0; m <= J; m++)
			{
				y[n] += b[m] * x[J + n - m];
			}
		}

		/* フィルタリング結果の連結 */
		for (n = 0; n < L / 2; n++)
		{
			pcm1.sL[(offset_ex / 2) + n] = y[2 * n];
			pcm1.sR[(offset_ex / 2) + n] = y[(2 * n) + 1];
		}
	}

	stereo_wave_write(&pcm1, "test2.WAV"); /* WAVEファイルにステレオの音データを出力する */

	stereo_wave_read(&pcm0, "test2.WAV"); /* WAVEファイルからステレオの音データを入力する */

	Hanning_window(w, N); /* ハニング窓 */

	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = N / 2 * frame;

		/* x(n)のFFT */
		for (n = 0; n < N/2; n++)
		{
			x_real[2 * n] = pcm0.sL[offset / 2 + n] * w[n];
			x_real[2 * n + 1] = pcm0.sR[offset / 2 + n] * w[n];
			x_imag[2 * n] = 0.0;
			x_imag[2 * n + 1] = 0.0;
		}
		FFT(x_real, x_imag, N);

		/* 振幅スペクトルと位相スペクトル */
		for (k = 0; k < N; k++)
		{
			A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
			if (x_imag[k] != 0.0 && x_real[k] != 0.0)
			{
				T[k] = atan2(x_imag[k], x_real[k]);
			}
		}

		/* スペクトルサブトラクション */
		for (k = 0; k < N; k++)
		{
			A[k] -= threshold;
			if (A[k] < 0.0)
			{
				A[k] = 0.0;
			}
		}

		for (k = 0; k < N; k++)
		{
			y_real[k] = A[k] * cos(T[k]);
			y_imag[k] = A[k] * sin(T[k]);
		}
		IFFT(y_real, y_imag, N);

		/* 加工結果の連結 */
		for (n = 0; n < N/2; n++)
		{
			pcm1.sL[offset / 2 + n] += y_real[2 * n];
			pcm1.sR[offset / 2 + n] += y_real[2 * n + 1];
		}
	}

	stereo_wave_write(&pcm1, "try.WAV"); /* WAVEファイルにステレオの音データを出力する */

	free(pcm0.sL); /* メモリの解放 */
	free(pcm0.sR); /* メモリの解放 */
	free(pcm1.sL); /* メモリの解放 */
	free(pcm1.sR); /* メモリの解放 */
	free(w); /* メモリの解放 */
	free(x_real); /* メモリの解放 */
	free(x_imag); /* メモリの解放 */
	free(A); /* メモリの解放 */
	free(T); /* メモリの解放 */
	free(y_real); /* メモリの解放 */
	free(y_imag); /* メモリの解放 */

	return 0;
}


