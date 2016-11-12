#pragma warning(disable:4996)

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "initfile.h"
#include "wave.h"
#include "window_function.h"
#include "sinc.h"
#include "fir_filter.h"
#include "fft.h"
#include "countSound.h"
#include "make_datfile.h"
#include "dfr_filter.h"
#include "spectrum.h"
#define OUT_FILE_NUM 4

int main(void)
{
	STEREO_PCM pcm0, bflap_pcm0, pcm1[OUT_FILE_NUM];
	IN_FILE inputfile;
	OUT_FILE outputfiles[OUT_FILE_NUM];
	int n, m, J, L, ife1, ife2, idelta, offset, frame, number_of_frame;
	double fe1, fe2, delta, *x, *xlp, *b, *wJ, *wN, *x_real, *x_imag, *x_lpre, *x_lpim, *y_real[OUT_FILE_NUM], *y_imag[OUT_FILE_NUM], *n_real, *n_imag;
	char *infilename, *outfilenames[OUT_FILE_NUM];


	//入力ファイルのデータ
	char file[] = "short1";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//BPフィルタのデータ
	ife1 = 1000;	/* エッジ周波数1 */
	ife2 = 2500;	/* エッジ周波数2 */
	idelta = 900;	/* 遷移帯域幅 */
					//フレームの長さ
	L = 1024;
	//ノイズのデータ
	int str_frame = 20;
	int fin_frame = 100;


	//ファイルの初期化処理
	init_infiles(&inputfile, directory, file, filetype);
	infilename = (char*)calloc(getinfsize(&inputfile), sizeof(char));
	getinfilename(&inputfile, infilename);
	init_outfiles(&inputfile, outputfiles, ife1, ife2, idelta);
	for (int i = 0; i < OUT_FILE_NUM; i++) {
		outfilenames[i] = (char*)calloc(getoutfsize(&outputfiles[i]), sizeof(char));
		getoutfilename(&outputfiles[i], outfilenames[i]);
	}

	/* WAVEファイルからステレオの音データを入力する */
	stereo_wave_read(&pcm0, infilename);

	/*各pcmの初期化*/
	spcmcpy(&bflap_pcm0, &pcm0);
	for (int i = 0; i < OUT_FILE_NUM; i++) {
		spcmcpy(&pcm1[i], &pcm0);
	}

	//BPフィルタｂを作る
	fe1 = (double)ife1 / (2 * pcm0.fs);
	fe2 = (double)ife2 / (2 * pcm0.fs);
	delta = (double)idelta / (2 * pcm0.fs);
	J = (int)(3.1 / delta + 0.5) - 1; /* 遅延器の数 */
	if (J % 2 == 1)
	{
		J++; /* J+1が奇数になるように調整する */
	}
	b = (double*)calloc((J + 1), sizeof(double));
	wJ = (double*)calloc((J + 1), sizeof(double));
	Hanning_window(wJ, (J + 1));
	FIR_BPF(fe1, fe2, J, b, wJ);

	/* フレームの数 */
	number_of_frame = (pcm0.length * 2 - L / 2) / (L / 2);

	/* メモリの確保 */
	x = (double*)calloc(L + J, sizeof(double));
	xlp = (double*)calloc(L + J, sizeof(double));
	x_real = (double*)calloc(L, sizeof(double));
	x_imag = (double*)calloc(L, sizeof(double));
	x_lpre = (double*)calloc(L, sizeof(double));
	x_lpim = (double*)calloc(L, sizeof(double));
	n_real = (double*)calloc(L, sizeof(double));
	n_imag = (double*)calloc(L, sizeof(double));
	wN = (double*)calloc(L, sizeof(double));
	Hanning_window(wN, L);
	/*出力フレームの初期化*/
	for (int i = 0; i < OUT_FILE_NUM; i++)
	{
		y_real[i] = (double*)calloc(L, sizeof(double));
		y_imag[i] = (double*)calloc(L, sizeof(double));
	}

	//推定ノイズの初期化
	for (n = 0; n < L; n++)
	{
		n_real[n] = 0.0;
		n_imag[n] = 0.0;
	}


	//フレーム単位の周波数領域での処理
	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = (L / 2) * frame;

		/*前・エッジ強調処理*/
		laplacian(x, xlp, L + J);

		/*BPF実行*/
		for (n = 0; n < (L + J) / 2; n++)/* 直前のフレームの後半のJサンプルをつけ加える */
		{
			if ((offset - J) / 2 + n < 0)
			{
				x[2 * n] = 0;
				x[(2 * n) + 1] = 0;
				xlp[2 * n] = 0;
				xlp[(2 * n) + 1] = 0;
			}
			else
			{
				x[2 * n] = pcm0.sL[(offset - J) / 2 + n];
				x[(2 * n) + 1] = pcm0.sR[(offset - J) / 2 + n];
				xlp[2 * n] = pcm0.sL[(offset - J) / 2 + n];
				xlp[(2 * n) + 1] = pcm0.sR[(offset - J) / 2 + n];
			}
		}
		for (n = 0; n < L; n++)/*初期化処理*/
		{
			x_real[n] = 0.0;
			x_lpre[n] = 0.0;
		}
		for (n = 0; n < L; n++)/* フィルタリング */
		{
			for (m = 0; m <= J; m++)
			{
				x_real[n] += b[m] * x[J + n - m];
				x_lpre[n] += b[m] * xlp[J + n - m];
			}
		}

		/*初期化処理*/
		for (n = 0; n < L / 2; n++)
		{
			x_imag[2 * n] = 0.0;
			x_imag[(2 * n) + 1] = 0.0;
			x_lpim[2 * n] = 0.0;
			x_lpim[(2 * n + 1)] = 0.0;
			for (int i = 0; i < OUT_FILE_NUM; i++) {
				y_real[i][n] = 0;
				y_imag[i][n] = 0;
			}
		}

		/*ノイズ設定*/
		if (frame > str_frame && frame < fin_frame)
		{
			setnoise(n_real, n_imag, x_real, L, frame);
			setnoise(n_real, n_imag, x_lpre, L, frame);
		}

		/* FFT */
		FFT(n_real, n_imag, L);
		FFT(x_real, x_imag, L);
		FFT(x_lpre, x_lpim, L);

		/*スペクトルサブトラクション(仮！！)*/
		do_subtruction(x_real, x_imag, n_real, n_imag, L);
		do_subtruction(x_lpre, x_lpim, n_real, n_imag, L);

		/* IFFT */
		IFFT(n_real, n_imag, L);
		IFFT(x_real, x_imag, L);
		IFFT(x_lpre, x_lpim, L);

		for (n = 0; n < L / 2; n++)
		{
			y_real[0][n] = x_real[n];
			y_real[2][n] = x_lpre[n];
		}
		/*後・エッジ強調処理*/
		laplacian(y_real[0], y_real[1], L);
		laplacian(y_real[2], y_real[3], L);

		/* フレームの連結 */
		for (n = 0; n < L / 2; n++)
		{
			if ((offset / 2) + n < pcm0.length)
			{
				pcm1[0].sL[(offset / 2) + n] += y_real[0][2 * n];
				pcm1[0].sR[(offset / 2) + n] += y_real[0][(2 * n) + 1];

				pcm1[2].sL[(offset / 2) + n] += y_real[2][2 * n];
				pcm1[2].sR[(offset / 2) + n] += y_real[2][(2 * n) + 1];
			}
		}
	}


	/* WAVEファイルにモノラルの音データを出力する */
	for (int i = 0; i < OUT_FILE_NUM; i++) {
		stereo_wave_write(&pcm1[i], outfilenames[i]);
	}

	printf("sec = %lf\n", (double)(number_of_frame * (L + 1)) / (double)(4 * pcm0.fs));

	/* メモリの解放 */
	free(infilename);
	free(pcm0.sL);
	free(pcm0.sR);
	free(bflap_pcm0.sL);
	free(bflap_pcm0.sR);
	free(b);
	free(wJ);
	free(wN);

	free(x);
	free(xlp);
	free(n_real);
	free(n_imag);
	free(x_real);
	free(x_imag);
	free(x_lpre);
	free(x_lpim);
	for (int i = 0; i < OUT_FILE_NUM; i++)
	{
		free(outfilenames[i]);
		free(pcm1[i].sL);
		free(pcm1[i].sR);
		free(y_real[i]);
		free(y_imag[i]);
	}

	printf("finish");

	return 0;
}
