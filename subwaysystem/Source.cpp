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

double getsec(int frame, int datanum, double fs);

int main(void)
{
	STEREO_PCM pcm0, pcm1;
	int m, n, J, k, N, L, offset, frame, number_of_frame, count, maxindex;
	double fe1, fe2, delta, *b, *b_real, *b_imag, *w, *x_real, *x_imag, *A, *y_real, *y_imag, max, calculate;

	stereo_wave_read(&pcm0, "short1.WAV"); /* WAVEファイルからステレオの音データを入力する */

	pcm1.fs = pcm0.fs; /* 標本化周波数 */
	pcm1.bits = pcm0.bits; /* 量子化精度 */
	pcm1.length = pcm0.length; /* 音データの長さ */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */

	fe1 = 1000.0 / (pcm0.fs * 2); /* エッジ周波数 */
	fe2 = 2000.0 / (pcm0.fs * 2); /* エッジ周波数 */
	delta = 800.0 / (pcm0.fs * 2); /* 遷移帯域幅 */

	J = (int)(3.1 / delta + 0.5) - 1; /* 遅延器の数 */
	if (J % 2 == 1)
	{
		J++; /* J+1が奇数になるように調整する */
	}

	L = 256; /* フレームの長さ */
	N = 512; /* DFTのサイズ */

	number_of_frame = 2 * pcm0.length / L; /* フレームの数 */
	/*printf("フレーム数：%d\n",number_of_frame);*/

	b = (double*)calloc((J + 1), sizeof(double)); /* メモリの確保 */
	w = (double*)calloc((J + 1), sizeof(double)); /* メモリの確保 */

	Hanning_window(w, (J + 1)); /* ハニング窓 */

	FIR_BPF(fe1, fe2, J, b, w); /* FIRフィルタの設計 */

	b_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	b_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	A = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */

	count = 0;
	max = A[0];
	calculate = pow(10, 2.25);

	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;
		max = 0.0;
		maxindex = 0;

		/* x(n)のFFT */
		for (n = 0; n < N; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
		}
		for (n = 0; n < L / 2; n++)
		{
			x_real[2 * n] = pcm0.sL[(offset / 2) + n];
			x_real[2 * n + 1] = pcm0.sR[(offset / 2) + n];
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

		/* 振幅スペクトル */
		for (k = 0; k < N; k++)
		{
			A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
		}

		/* フィルタリング */
		for (k = 0; k < N; k++)
		{
			y_real[k] = x_real[k] * b_real[k] - x_imag[k] * b_imag[k];
			y_imag[k] = x_imag[k] * b_real[k] + x_real[k] * b_imag[k];
		}

		IFFT(y_real, y_imag, N);

		/* フィルタリング結果の連結 */
		for (n = 0; n < L; n++)
		{
			if (offset / 2 + n < pcm1.length)
			{
				pcm1.sL[offset / 2 + n] += y_real[2 * n];
				pcm1.sR[offset / 2 + n] += y_real[2 * n + 1];
			}
		}

		for (k = 0; k < N; k++)
		/*{
			if (A[k] > max)
			{
				max = A[k];
				maxindex = k;
			}
			printf("A[0]の値:%lf\n",A[0]);
		}*/
		if (A[k] > calculate*A[0])
		{
				count++;
				/*printf("frame = %d, time = %lf, max:A[%d]=%lf\n",frame,getsec(frame,L,pcm0.fs), maxindex, max);*/
		}

	}
	printf("回数：%d", count);
	/*printf("fin_time = %lf\n",getsec(number_of_frame,L,pcm0.fs));*/

	stereo_wave_write(&pcm1, "try.WAV"); /* WAVEファイルにステレオの音データを出力する */

	free(pcm0.sL); /* メモリの解放 */
	free(pcm0.sR); /* メモリの解放 */
	free(pcm1.sL); /* メモリの解放 */
	free(pcm1.sR); /* メモリの解放 */
	free(w); /* メモリの解放 */
	free(x_real); /* メモリの解放 */
	free(x_imag); /* メモリの解放 */
	free(A); /* メモリの解放 */
	free(y_real); /* メモリの解放 */
	free(y_imag); /* メモリの解放 */

	return 0;
}


double getsec(int frame, int datanum, double fs) {
	double sec;

	sec = ((double)(datanum*(frame + 1)) / (double)(4 * fs));

	return 0;
}