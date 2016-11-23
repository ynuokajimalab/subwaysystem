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
	STEREO_PCM pcm0, cut_pcm0, noise_pcm0, pcm1[OUT_FILE_NUM];
	IN_FILE inputfile[2];
	OUT_FILE outputfiles[OUT_FILE_NUM];
	int n, m, k, Jb, N, L, offset, frame, number_of_frame, filternum;
	int str_index, fin_index, ife1, ife2, idelta;
	int noise_index1, noise_index2, noise_size;
	int index, maxi, mini, histgramsize;
	double max, min, data, dy, counttime, s;
	int maxindex, minindex, secondindex, bflag;
	int cflag, count, sleepframe;
	long Ip, Jp;
	double fe1, fe2, delta;
	double str_time, fin_time;
	double noise_time1, noise_time2, noiseaverage, noiseSD_t, noiseSD_f, alpha;
	double dn, hd, *histgram1, *histgram2;
	double threshold, sleeptime;
	double *b_real, *b_imag, *bb, *wJb, *wN, *x_real, *x_imag, *x_lpre, *x_lpim, *y_real[OUT_FILE_NUM], *y_imag[OUT_FILE_NUM];
	double *n_real, *n_imag, *x_amp, *x_phas, *n_amp, *n_phas;
	double fc, Q, g, a[3], bp[3], **A, **B, *Ael, *Bel;
	char *infilename[2], *outfilenames[OUT_FILE_NUM], *datafile;


	//入力ファイルのデータ
	char file[] = "short2";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//フレームの長さ
	N = 1024;
	L = 512;  //実際のデータサイズ
	//ノイズのデータ
	noise_time1 = 4.0;
	noise_time2 = 8.0;
	alpha = 5.2;
	dn = 0.005;
	s = 0.04;
	//
	bflag = 0;
	cflag = 0;
	sleeptime = 0.5;
	counttime = 0.0;

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
	printf("ファイルを読み込みました\n");
	printf("総データ数：%d[個]\n", (pcm0.length) * 2);


	//サイズNの窓関数作成
	wN = (double*)calloc(N, sizeof(double));
	Hanning_window(wN, N);


	//ノイズ設定
	noise_index1 = getstereoindex(noise_time1, pcm0.fs);
	noise_index2 = getstereoindex(noise_time2, pcm0.fs);
	noise_size = noise_index2 - noise_index1;
	printf("\n----ノイズに関するデータ----\n");
	printf("ノイズ開始時刻：%lf, ノイズ終了時刻：%lf\n", noise_time1, noise_time2);
	printf("ノイズ開始データ = pcm[%d], ノイズ終了データ = pcm[%d]\n", noise_index1, noise_index2);
	setspcm(&noise_pcm0, &pcm0, noise_size / 2);
	getpcm(&pcm0, &noise_pcm0, noise_index1, noise_index2);

	n_real = (double*)calloc(noise_size, sizeof(double));
	n_imag = (double*)calloc(noise_size, sizeof(double));
	n_amp = (double*)calloc(noise_size, sizeof(double));
	n_phas = (double*)calloc(noise_size, sizeof(double));

	//推定ノイズの初期化
	for (n = 0; n < noise_size; n++)
	{
		n_amp[n] = 0.0;
		n_phas[n] = 0.0;
		n_real[n] = 0.0;
		n_imag[n] = 0.0;
	}
	setnoise(n_real, n_imag, n_amp, n_phas, noise_pcm0);
	//noiseaverage = getsoundmyu(n_real,noise_size);
	noiseSD_t = getSD(n_real, 0, noise_size);
	threshold = noiseSD_t * alpha;
	printf("平均：%lf\n", noiseaverage);
	printf("標準偏差：%lf\n", noiseSD_t);
	printf("閾値：%lf\n", threshold);

	/*各pcmの初期化*/
	for (int i = 0; i < OUT_FILE_NUM; i++) {
		spcmcpy(&pcm1[i], &pcm0);
	}

	/* フレームの数 */
	number_of_frame = (pcm0.length * 2 - L / 2) / (L / 2);
	printf("フレーム数：%d\n", number_of_frame);

	printf("\n----カウントに関するデータ----\n");
	sleepframe = getframe(sleeptime, L, pcm0.fs);
	count = 0;
	printf("閾値：%lf\nカウント休止時間：%lf[s]\nカウント休止フレーム数：%d\n", threshold, sleeptime, sleepframe);


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
		offset = (L / 2) * frame;

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
			for (int i = 0; i < OUT_FILE_NUM; i++) {
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


		for (int i = 0; i < OUT_FILE_NUM; i++) {
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

		//音のカウント
		if (bflag == 1) {
			secondindex = getMaxindex(y_real[0], N);
			dy = delta - y_real[0][secondindex];
			if (dy > 0.15) {
				double temptime;
				count++;
				temptime = getsecond(frame, L, pcm0.fs);
				if (sleeptime > (temptime - counttime)) {
					sleeptime -= (temptime - counttime) * s;
				}
				else if (sleeptime < (temptime - counttime)) {
					sleeptime += (temptime - counttime) * s;
				}
				counttime = temptime;
				sleepframe = getframe(sleeptime, L, pcm0.fs);
				printf("count = %d, frame = %d, time = %lf\n", count, frame, counttime);
				printf("カウント休止時間：%lf[s]\nカウント休止フレーム数：%d\n", sleeptime, sleepframe);
				printf("max:sound[%d] = %lf, min:sound[%d] = %lf\n\n", maxindex, max, minindex, min);
				cflag = sleepframe;

			}
			bflag = 0;
		}
		if (cflag == 0) {


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

			if ((delta > threshold)) {
				//printf("%d：閾値を超えました\n", frame);
				bflag = 1;
			}
		}
		else
		{
			cflag--;
		}

		/* フレームの連結 */
		for (n = 0; n < N / 2; n++)
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

	printf("解析時間： %lf\n", getsecond(number_of_frame, L, pcm0.fs));
	printf("----解析終了----\n");


	/* WAVEファイルにモノラルの音データを出力する */
	for (int i = 0; i < OUT_FILE_NUM; i++) {
		stereo_wave_write(&pcm1[i], outfilenames[i]);
	}
	stereo_wave_write(&cut_pcm0, infilename[1]);
	printf("ファイルを出力しました\n");


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
	for (int i = OUT_FILE_NUM - 1; i >= 0; i--)
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
