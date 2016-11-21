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
#include "./headerfiles/iir_filter.h"
#include "./headerfiles/fft.h"
#include "./headerfiles/countSound.h"
#include "./headerfiles/make_datfile.h"
#include "./headerfiles/dfr_filter.h"
#include "./headerfiles/spectrum.h"
#include "./headerfiles/frame-time.h"
#include "./headerfiles/statistics.h"
#define OUT_FILE_NUM 4

int main(void)
{
	STEREO_PCM pcm0, cut_pcm0,noise_pcm0, pcm1[OUT_FILE_NUM];
	IN_FILE inputfile[2];
	OUT_FILE outputfiles[OUT_FILE_NUM];
	SPECTROGRAM **sp, *spel;
	int n, m, k, Jb, N, L, offset, frame, number_of_frame, filternum;
	int str_index, fin_index,ife1, ife2, idelta;
	int noise_index1, noise_index2, noise_size;
	int cflag, count, sleepframe;
	long Ip, Jp;
	double fe1, fe2, delta;
	double str_time, fin_time;
	double noise_time1, noise_time2,noiseaverage, noiseSD_t,noiseSD_f,alpha;
	double threshold, sleeptime;
	double *b_real, *b_imag, *bb, *wJb, *wN, *x_real, *x_imag, *x_lpre, *x_lpim, *y_real[OUT_FILE_NUM], *y_imag[OUT_FILE_NUM];
	double *n_real, *n_imag, *x_amp, *x_phas, *n_amp, *n_phas;
	double fc, Q, g, a[3], bp[3], **A, **B, *Ael, *Bel;
	char *infilename[2], *outfilenames[OUT_FILE_NUM],*datafile;


	//入力ファイルのデータ
	char file[] = "short1";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//BPフィルタのデータ
	ife1 = 1000;	/* エッジ周波数1 */
	ife2 = 2000;	/* エッジ周波数2 */
	idelta = 800;	/* 遷移帯域幅 */
	//ピーキングフィルタのデータ
	Ip = 2; /* 遅延器の数 */
	Jp = 2; /* 遅延器の数 */
	Q = 2.0; /* クオリティファクタ（一定） */	
	filternum = 1024;
	//フレームの長さ
	N = 1024; 
	L = 512;  //実際のデータサイズ
	//解析する時間区間
	str_time = 0.0;
	fin_time = 7.0;
	//ノイズのデータ
	noise_time1 = 1.0;
	noise_time2 = 7.3;
	alpha = 3.0;
	//
	cflag = 0;
	sleeptime = 0.5;

	printf("ファイル名：%s\n", file);
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


	//ノイズ設定
	noise_index1 = getstereoindex(noise_time1,pcm0.fs);
	noise_index2 = getstereoindex(noise_time2, pcm0.fs);
	noise_size = noise_index2 - noise_index1;
	printf("\n----ノイズに関するデータ----\n");
	printf("ノイズ開始時刻：%lf, ノイズ終了時刻：%lf\n", noise_time1, noise_time2);
	printf("ノイズ開始データ = pcm[%d], ノイズ終了データ = pcm[%d]\n",noise_index1,noise_index2);
	setspcm(&noise_pcm0, &pcm0, noise_size / 2);
	getpcm(&pcm0, &noise_pcm0, noise_index1, noise_index2);

	n_real = (double*)calloc(2*noise_pcm0.length, sizeof(double));
	n_imag = (double*)calloc(2 * noise_pcm0.length, sizeof(double));
	n_amp = (double*)calloc(2 * noise_pcm0.length, sizeof(double));
	n_phas = (double*)calloc(2 * noise_pcm0.length, sizeof(double));

	//推定ノイズの初期化
	for (n = 0; n < noise_size; n++)
	{
		n_amp[n] = 0.0;
		n_phas[n] = 0.0;
		n_real[n] = 0.0;
		n_imag[n] = 0.0;
	}
	setnoise(n_real, n_imag, n_amp, n_phas, noise_pcm0);
	noiseaverage = getsoundmyu(n_real,noise_size);
	noiseSD_t = getSD(n_real,noiseaverage,noise_size);
	threshold = noiseaverage + noiseSD_t * alpha;
	printf("平均：%lf\n", noiseaverage);
	printf("標準偏差：%lf\n", noiseSD_t);
	printf("閾値：%lf\n", threshold);

	//解析区間の設定
	str_index = getstereoindex(str_time,pcm0.fs);
	fin_index = getstereoindex(fin_time, pcm0.fs);
	printf("\n----解析データ----\n");
	printf("開始時刻：%lf 終了時刻: %lf\n", str_time,fin_time);
	printf("開始データ：pcm[%d], 終了データ：pcm[%d]\n", str_index, fin_index);
	setspcm(&cut_pcm0, &pcm0, (fin_index - str_index) / 2);
	getpcm(&pcm0,&cut_pcm0,str_index,fin_index);

	/*各pcmの初期化*/

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
	Jb = (int)(3.1 / delta + 0.5) - 1; /* 遅延器の数 */
	if (Jb % 2 == 1)
	{
		Jb++; /* Jb+1が奇数になるように調整する */
	}
	bb = (double*)calloc((Jb + 1), sizeof(double));
	wJb = (double*)calloc((Jb + 1), sizeof(double));
	Hanning_window(wJb, (Jb + 1));
	FIR_BPF(fe1, fe2, Jb, bb, wJb);

	//ピーキングフィルタを作る
	A = (double**)calloc(N ,sizeof(double*));
	B = (double**)calloc(N, sizeof(double*));
	Ael = (double*)calloc(N * 3, sizeof(double));
	Bel = (double*)calloc(N * 3, sizeof(double));
	for (k = 0; k<N; k++) {
		A[k] = Ael + k * Ip;
		B[k] = Bel + k * Jp;
	}
	for (k = 0; k < filternum; k++) {
		fc = k / pcm0.fs; /* 中心周波数 */
		g = -n_amp[k];
		/*int index = getMaxindex(n_amp, N);
		g = n_amp[index];*/
		IIR_peaking(fc, Q, g, a, bp); /* IIRフィルタの設計 */
		for (m = 0; m <= Ip; m++)
		{
			A[k][m] = a[m];
		}
		for (m = 0; m <= Jp; m++)
		{
			B[k][m] = bp[m];
		}
	}


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

	////データファイルの作成
	//datafile = (char*)calloc(10,sizeof(char));
	//strcpy(datafile,file);
	//strcat(datafile,".dat");
	//make_datfile(datafile);
	//sp = (SPECTROGRAM**)calloc(number_of_frame, sizeof(SPECTROGRAM*));
	//spel = (SPECTROGRAM*)calloc(N * number_of_frame, sizeof(SPECTROGRAM));
	//for (frame = 0; frame < number_of_frame; frame++) {
	//	sp[frame] = spel + frame * N;
	//}


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

		for (m = 0; m <= Jb; m++)
		{
			b_real[m] = bb[m];
		}

		/*前・エッジ強調処理*/
		laplacian(x_real, x_lpre, N);

		for (int i = 0; i < N; i++) {//窓関数をかける
			x_real[i] *= wN[i];
			x_lpre[i] *= wN[i];
		}
		FFT(x_real, x_imag, N);
		FFT(x_lpre, x_lpim, N);

		///*スペクトルサブトラクション(仮！！)*/
		//for (int k = 0; k < L; k++) {//初期化
		//	x_amp[k] = 0.0;
		//	x_phas[k] = 0.0;
		//}
		//getampphase(x_real, x_imag, x_amp, x_phas, N);
		//subtruction(x_amp, n_amp, N);
		//getrealimage(x_real, x_imag, x_amp, x_phas, N);
		//for (int k = 0; k < L; k++) {//初期化
		//	x_amp[k] = 0.0;
		//	x_phas[k] = 0.0;
		//}
		//getampphase(x_lpre, x_lpim, x_amp, x_phas, N);
		//subtruction(x_amp, n_amp, N);
		//getrealimage(x_lpre, x_lpim, x_amp, x_phas, N);
		///*for (int k = 0; k < L; k++) {
		//x_real[k] *= wN[k];
		//x_imag[k] *= wN[k];
		//}*/

		//データをセーブ
		for (n = 0; n < N; n++)
		{
			for (int i = 0; i < OUT_FILE_NUM;i++) {
			y_real[i][n] = x_real[n];	
			y_imag[i][n] = x_imag[n];
			}
		}

	/*    BPF 実行
		FFT(b_real, b_imag, N);
		for (n = 0; n < N; n++) {
			y_real[0][n] = x_real[n] * b_real[n] - x_imag[n] * b_imag[n];
			y_imag[0][n] = x_imag[n] * b_real[n] + x_real[n] * b_imag[n];
			
			y_real[2][n] = x_lpre[n] * b_real[n] - x_lpim[n] * b_imag[n];
			y_imag[2][n] = x_lpim[n] * b_real[n] + x_lpre[n] * b_imag[n];
		}
*/


		for (int i = 0; i < OUT_FILE_NUM;i++) {
			IFFT(y_real[i], y_imag[i], N);
		}

		//データを更新
		for (n = 0; n < N; n++)
		{
			x_real[n] = y_real[0][n];
			//x_imag[n] = y_imag[0][n];
			x_lpre[n] = y_real[2][n];
			//x_lpim[n] = y_imag[2][n];
		}

		////イコライザ
		//for (k = 0; k < filternum; k++)	{
		//	for (n = 0; n < N; n++){
		//		y_real[0][n] = 0.0;
		//	}
		//	for (n = 0; n < N; n++){
		//		for (m = 0; m <= Jp; m++){
		//			if (n - m >= 0){
		//				y_real[0][n] += B[k][m] * x_real[n - m];
		//			}
		//		}
		//		for (m = 1; m <= Ip; m++)
		//		{
		//			if (n - m >= 0)
		//			{
		//				y_real[0][n] += -A[k][m] * x_real[n - m];
		//			}
		//		}
		//	}
		//}

		//データを更新
		for (n = 0; n < N; n++)
		{
			x_real[n] = y_real[0][n];
			//x_imag[n] = y_imag[0][n];
			x_lpre[n] = y_real[2][n];
			//x_lpim[n] = y_imag[2][n];
		}


		/*後・エッジ強調処理*/
		laplacian(y_real[0], y_real[1], N);
		laplacian(y_real[2], y_real[3], N);

		//音のカウント
		if (cflag == 0) {

			int maxindex, minindex,index;
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
			//delta = max - min;
			if (max > -min) {
				delta = max;
				index = maxindex;
			}
			else {
				delta = min;
				index = minindex;
			}


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

		////.datファイル作成
		//double time = getsecond(frame,L,pcm0.fs);
		//double finterval = (pcm0.fs / N);
		//double db;
		//for (n = 0; n < N; n++) {
		//	sp[frame][n].time = time;
		//	sp[frame][n].freequency = (int)(finterval * n);
		//	//printf("x_amp[%d] = %lf\n",n,x_amp[n]);
		//	sp[frame][n].power =x_amp[n];
		//}
		//postscript_datfile(datafile, sp[frame], N);

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
	//free(spel);
	//free(sp);
	free(x_phas);
	free(x_amp);
	free(x_lpim);
	free(x_lpre);
	free(x_imag);
	free(x_real);
	free(b_imag);
	free(b_real);
	free(Bel);
	free(Ael);
	free(B);
	free(A);
	free(wJb);
	free(bb);
	free(n_phas);
	free(n_amp);
	free(n_imag);
	free(n_real);
	free(wN);
	for (int i = OUT_FILE_NUM-1; i >= 0; i--)
	{
		free(outfilenames[i]);
		free(pcm1[i].sL);
		free(pcm1[i].sR);
		free(y_real[i]);
		free(y_imag[i]);
	}
	free(infilename[1]);
	free(infilename[0]);
	free(pcm0.sR);
	free(pcm0.sL);
	free(cut_pcm0.sR);
	free(cut_pcm0.sL);

	printf("finish\n");

	return 0;
}
