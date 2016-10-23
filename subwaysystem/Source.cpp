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
#include "make_datfile.h"


int main(void)
{
	STEREO_PCM pcm0, pcm1;
	int n, m, k, J, L, N, offset, frame, number_of_frame,countsound,compressionrate,filehasdata;
	double fe1,fe2, delta, threshold, *b, *w, *b_real, *b_imag, *x_real, *x_imag, *y_real, *y_imag;


	stereo_wave_read(&pcm0, "short1.wav"); /* WAVEファイルからステレオの音データを入力する */

	pcm1.fs = pcm0.fs; /* 標本化周波数 */
	pcm1.bits = pcm0.bits; /* 量子化精度 */
	pcm1.length = pcm0.length; /* 音データの長さ */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* メモリの確保 */

	fe1 = 150.0 / pcm0.fs; /* エッジ周波数1 */
	fe2 = 450.0 / pcm0.fs; /* エッジ周波数2 */
	delta = 450.0 / pcm0.fs; /* 遷移帯域幅 */

	J = (int)(3.1 / delta + 0.5) - 1; /* 遅延器の数 */
	if (J % 2 == 1)
	{
		J++; /* J+1が奇数になるように調整する */
	}

	b = (double*) calloc((J + 1), sizeof(double)); /* メモリの確保 */
	w = (double*)calloc((J + 1), sizeof(double)); /* メモリの確保 */

	Hanning_window(w, (J + 1)); /* ハニング窓 */

	FIR_BPF(fe1,fe2, J, b, w); /* FIRフィルタの設計 */

	L = 256; /* フレームの長さ */
	N = 512; /* DFTのサイズ */

	countsound = 0; /* 音が閾値を超えた回数*/
	threshold = 10.0;/* 閾値*/

	filehasdata = 0; /*ファイルを新しく作るか決定するフラグ*/
	compressionrate = 500; /*データの圧縮率(-倍)*/


	number_of_frame = pcm0.length*2 / L; /* フレームの数 */


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
		for (n = 0; n < L/2; n++)
		{
			x_real[2 * n] = pcm0.sL[(offset/2) + n];
			x_real[(2 * n) + 1] = pcm0.sR[(offset / 2) + n];
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

		if (filehasdata != 0) {
			postscript_datfile("allsound.dat", y_real, N,offset);
			postscript_datfile("orgsound.dat", x_real, N, offset);

		}
		else {
			make_datfile("allsound.dat",y_real,N);
			make_datfile("orgsound.dat", x_real, N);
			filehasdata = 1;
		}

		if (frame == number_of_frame/2)
			make_datfile("onepointsound.dat",y_real,N);


		/* フィルタリング結果の連結 */
		for (n = 0; n < L * 2/2; n++)
		{
			if ((offset/2) + n < pcm1.length)
			{	
				pcm1.sL[(offset/2) + n] += y_real[2 * n];			
				pcm1.sR[(offset/2) + n] += y_real[(2*n)+1];
			}
		}
	}
	double max1 = sampling_data("allsound.dat", "samplingsound.dat", compressionrate);
	normalize_data("samplingsound.dat","nomalizesound.dat",max1);

	stereo_wave_write(&pcm1, "out_short1_BPFver.wav"); /* WAVEファイルにモノラルの音データを出力する */

	free(pcm0.sL); /* メモリの解放 */
	free(pcm0.sR); /* メモリの解放 */
	free(pcm1.sL); /* メモリの解放 */
	free(pcm1.sR); /* メモリの解放 */
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
