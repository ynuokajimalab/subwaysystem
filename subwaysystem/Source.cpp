#pragma warning(disable:4996)

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "./headerfiles/initfile.h"
#include "./headerfiles/wave.h"
#include "./headerfiles/window_function.h"
#include "./headerfiles/sinc.h"
#include "./headerfiles/fir_filter.h"
#include "./headerfiles/fft.h"
#include "./headerfiles/countSound.h"
#include "./headerfiles/make_datfile.h"
#include "./headerfiles/dfr_filter.h"
#include "./headerfiles/spectrum.h"
#include "./headerfiles/frame-time.h"
#define OUT_FILE_NUM 4

int main(void)
{
	STEREO_PCM pcm0, cut_pcm0, pcm1[OUT_FILE_NUM];
	IN_FILE inputfile[2];
	OUT_FILE outputfiles[OUT_FILE_NUM];
	int n, m, J,N, L,str_index,fin_index, ife1, ife2, idelta, offset, frame, number_of_frame, noise_frame1, noise_frame2, cflag,count, sleepframe;
	double fe1, fe2, delta,str_time,fin_time,noise_time1, noise_time2, threshold, sleeptime,*b_real,*b_imag, *b, *wJ,*wN, *x_real, *x_imag, *x_lpre, *x_lpim, *y_real[OUT_FILE_NUM], *y_imag[OUT_FILE_NUM], *n_real, *n_imag, *x_amp, *x_phas, *n_amp, *n_phas;
	char *infilename[2], *outfilenames[OUT_FILE_NUM];


	//入力ファイルのデータ
	char file[] = "short2";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//BPフィルタのデータ
	ife1 = 0000;	/* エッジ周波数1 */
	ife2 = 2000;	/* エッジ周波数2 */
	idelta = 600;	/* 遷移帯域幅 */
	//フレームの長さ
	N = 1024; 
	L = 512;  //実際のデータサイズ
	//解析する時間区間
	str_time = 0.0;
	fin_time = 12.5;
	//ノイズのデータ
	noise_time1 = 4.4;
	noise_time2 = 4.7;
	//閾値
	threshold = 0.2;
	cflag = 0;
	sleeptime = 0.5;

	printf("バンドパス： %d-%d[Hz], 遷移帯域幅：%d[Hz]\n", ife1, ife2, idelta);

	
	//ファイルの初期化処理
	init_infiles(inputfile, directory, file, filetype);
	infilename[0] = (char*)calloc(getinfsize(&inputfile[0]), sizeof(char));
	infilename[1] = (char*)calloc(4 + getinfsize(&inputfile[1]), sizeof(char));
	getinfilename(inputfile, infilename);
	init_outfiles(&inputfile[0], outputfiles, ife1, ife2, idelta);
	for (int i = 0; i < OUT_FILE_NUM; i++) {
		outfilenames[i] = (char*)calloc(getoutfsize(&outputfiles[i]), sizeof(char));
		getoutfilename(&outputfiles[i], outfilenames[i]);
	}

	/* WAVEファイルからステレオの音データを入力する */
	stereo_wave_read(&pcm0, infilename[0]);
	printf("総データ数：%d[個]\n", (pcm0.length) * 2);


	//サイズNの窓関数作成
	wN = (double*)calloc(N, sizeof(double));
	Hanning_window(wN, N);


	n_real = (double*)calloc(N, sizeof(double));
	n_imag = (double*)calloc(N, sizeof(double));
	n_amp = (double*)calloc(N, sizeof(double));
	n_phas = (double*)calloc(N, sizeof(double));

	//推定ノイズの初期化
	for (n = 0; n < N; n++)
	{
		n_amp[n] = 0.0;
		n_phas[n] = 0.0;
		n_real[n] = 0.0;
		n_imag[n] = 0.0;
	}
	//ノイズ設定
	noise_frame1 = getframe(noise_time1,L,pcm0.fs);
	noise_frame2 = getframe(noise_time2, L, pcm0.fs);
	printf("\n----ノイズに関するデータ----\n");
	printf("ノイズ開始時刻：%lf, ノイズ終了時刻：%lf\n", noise_time1, noise_time2);
	printf("ノイズ開始フレーム = %d, ノイズ終了フレーム = %d\n",noise_frame1,noise_frame2);
	setnoise(n_real, n_imag, n_amp, n_phas, wN,pcm0, N,L, noise_frame1, noise_frame2);

	//解析区間の設定
	str_index = getstereoindex(str_time,pcm0.fs);
	fin_index = getstereoindex(fin_time, pcm0.fs);
	printf("\n----解析データ----\n");
	printf("開始時間：%lf 終了時間 : %lf\n", str_time,fin_time);
	printf("開始データ：pcm[%d], 終了データ：pcm[%d]\n", str_index, fin_index);


	/*各pcmの初期化*/
	//spcmcpy(&cut_pcm0, &pcm0);
	cut_pcm0.fs = pcm0.fs; /* 標本化周波数 */
	cut_pcm0.bits = pcm0.bits; /* 量子化精度 */
	cut_pcm0.length = (fin_index - str_index)/2; /* 音データの長さ */
	cut_pcm0.sL = (double*)calloc(cut_pcm0.length, sizeof(double)); /* メモリの確保 */
	cut_pcm0.sR = (double*)calloc(cut_pcm0.length, sizeof(double)); /* メモリの確保 */
	getpcm(&pcm0,&cut_pcm0,str_time,fin_time);

	pcm0 = cut_pcm0;

	for (int i = 0; i < OUT_FILE_NUM; i++) {
		spcmcpy(&pcm1[i], &pcm0);
	}

	/* フレームの数 */
	number_of_frame = (pcm0.length * 2 - L/2) / (L/2);
	printf("フレーム数：%d\n", number_of_frame);

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


	sleepframe = getframe(sleeptime,L,pcm0.fs);
	count = 0;
	printf("\n----カウントに関するデータ----\n");
	printf("閾値：%lf\nカウント休止時間：%lf[s]\nカウント休止フレーム数：%d\n", threshold,sleeptime, sleepframe);
	

	/* メモリの確保 */
	b_real = (double*)calloc(N, sizeof(double));
	b_imag = (double*)calloc(N, sizeof(double));
	x_real = (double*)calloc(N, sizeof(double));
	x_imag = (double*)calloc(N, sizeof(double));
	x_lpre = (double*)calloc(N, sizeof(double));
	x_lpim = (double*)calloc(N, sizeof(double));
	x_amp = (double*)calloc(N, sizeof(double));
	x_phas = (double*)calloc(N, sizeof(double));
	/*出力フレームの初期化*/
	for (int i = 0; i < OUT_FILE_NUM; i++)
	{
		y_real[i] = (double*)calloc(N, sizeof(double));
		y_imag[i] = (double*)calloc(N, sizeof(double));
	}


	printf("\n----解析開始----\n");
	//フレーム単位の周波数領域での処理
	for (frame = 0; frame < number_of_frame; frame++)
	{
		//printf("---------------%dフレーム---------------------\n",frame);
		offset = (L/2) * frame;

		/*データを0で初期化*/
		for (n = 0; n < N; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
			x_lpre[n] = 0.0;
			x_lpim[n] = 0.0;
			b_real[n] = 0.0;
			b_imag[n] = 0.0;
			for (int i = 0; i < OUT_FILE_NUM; i++) {
				y_real[i][n] = 0;
				y_imag[i][n] = 0;
			}
		}

		//データ入力
		for (n = 0; n < L / 2; n++)
		{
			x_real[2 * n] = pcm0.sL[(offset / 2) + n];
			x_real[(2 * n) + 1] = pcm0.sR[(offset / 2) + n];
		}

		for (m = 0; m <= J; m++)
		{
			b_real[m] = b[m];
		}

		/*前・エッジ強調処理*/
		laplacian(x_real, x_lpre, N);

		for (int i = 0; i < N; i++) {//窓関数をかける
			x_real[i] *= wN[i];
			x_lpre[i] *= wN[i];
		}
		FFT(x_real, x_imag, N);
		FFT(x_lpre, x_lpim, N);

		/*スペクトルサブトラクション(仮！！)*/
		for (int k = 0; k < L; k++) {//初期化
			x_amp[k] = 0.0;
			x_phas[k] = 0.0;
		}
		getampphase(x_real, x_imag, x_amp, x_phas, N);
		subtruction(x_amp, n_amp, N);
		getrealimage(x_real, x_imag, x_amp, x_phas, N);
		/*for (int k = 0; k < L; k++) {
		x_real[k] *= wN[k];
		x_imag[k] *= wN[k];
		}*/

		//データをセーブ
		for (n = 0; n < N; n++)
		{
			y_real[1][n] = x_real[n];
			y_real[3][n] = x_lpre[n];
		}

		//BPF 実行
		FFT(b_real, b_imag, N);
		for (n = 0; n < N; n++) {
			y_real[0][n] = x_real[n] * b_real[n] - x_imag[n] * b_imag[n];
			y_imag[0][n] = x_imag[n] * b_real[n] + x_real[n] * b_imag[n];
			
			y_real[2][n] = x_lpre[n] * b_real[n] - x_lpim[n] * b_imag[n];
			y_imag[2][n] = x_lpim[n] * b_real[n] + x_lpre[n] * b_imag[n];


		}

		for (int i = 0; i < OUT_FILE_NUM;i++) {
			IFFT(y_real[i], y_imag[i], N);
		}

		///*後・エッジ強調処理*/
		//laplacian(y_real[0], y_real[1], N);
		//laplacian(y_real[2], y_real[3], N);

		//音のカウント
		if (cflag == 0) {

			int maxindex, minindex;
			double max, min, data, delta;
	
			max = -1.0;
			min = 1.0;
			for (int i = 0; i < N; i++) {
				data = y_real[0][i];
				if (max < data) {
					max = data;
					maxindex = i;
				}
				if (min > data) {
					min = data;
					minindex = i;
				}
			}
			delta = max - min;
			if (delta > threshold) {
				count ++;
				printf("count = %d, frame = %d, time = %lf\n", count,frame, getsecond(frame, L,pcm0.fs));
				printf("max:sound[%d] = %lf, min:sound[%d] = %lf\n\n", maxindex, max, minindex, min);
				cflag = sleepframe;
			}
		}else 
		{
			cflag--;
		}


		/* フレームの連結 */
		for (n = 0; n < N/2; n++)
		{
			if ((offset / 2) + n < pcm0.length)
			{
				for (int i = 0; i < OUT_FILE_NUM; i++)
				{
				pcm1[i].sL[(offset / 2) + n] += y_real[i][2 * n];
				pcm1[i].sR[(offset / 2) + n] += y_real[i][(2 * n) + 1];
				
				}

				pcm1[2].sL[(offset / 2) + n] += y_real[2][2 * n];
				pcm1[2].sR[(offset / 2) + n] += y_real[2][(2 * n) + 1];
			}
		}
	}

	printf("解析時間： %lf\n",getsecond(number_of_frame,L,pcm0.fs));
	printf("----解析終了----\n");


	/* WAVEファイルにモノラルの音データを出力する */
	for (int i = 0; i < OUT_FILE_NUM; i++) {
		stereo_wave_write(&pcm1[i], outfilenames[i]);
	}
	stereo_wave_write(&cut_pcm0,infilename[1]);

	/* メモリの解放 */
	free(infilename[0]);
	free(infilename[1]);
	free(pcm0.sL);
	free(pcm0.sR);
	//free(cut_pcm0.sL);
	//free(cut_pcm0.sR);
	free(b);
	free(wJ);
	free(wN);
	free(x_amp);
	free(x_phas);
	free(n_amp);
	free(n_phas);
	
	free(b_real);
	free(b_imag);
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

	printf("finish\n");

	return 0;
}
