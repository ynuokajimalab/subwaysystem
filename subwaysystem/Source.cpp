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

int main(void)
{
	STEREO_PCM org_pcm, noise_pcm, out_pcm;
	FILENAME inputfile,outputfile;
	size_t outfsize;
	int n, N, L, offset, frame, number_of_frame;
	int noise_index1, noise_index2, noise_size;
	int maxindex, minindex, representive_index,secondindex, bflag, cflag, count, sleepframe;
	double *x_real, *x_imag, *y_real,*y_imag;
	double noise_time1, noise_time2, noiseSD_t,alpha, *n_real, *n_imag;
	double threshold, sleeptime, max, min, data, dy,judgedy, counttime, s,temptime,representive_value;
	char *infilename, *outfilename;


	//入力ファイルのデータ
	char orgfile[] = "short3";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//フレームの長さ
	N = 1024;
	L = 512;  //実際のデータサイズ
	//ノイズのデータ
	noise_time1 = 0.0;
	noise_time2 = 7.0;
	sleeptime = 0.3;
	judgedy = 0.04;
	alpha = 4.9;
	s = 0.02;


	//データの初期化
	bflag = 0;	//前フレームで閾値を超えたか判定
	cflag = 0;	//カウントをチェックするかどうか判定
	counttime = 0.0;	//前のカウント時刻

	//ファイルの初期化処理
	inputfile.filename = (char*)calloc(20, sizeof(char));
	outputfile.filename = (char*)calloc(20, sizeof(char));
	init_file(&inputfile, directory, orgfile, filetype);
	init_file(&outputfile, directory, orgfile, filetype);
	infilename = (char*)calloc(getfnsize(inputfile), sizeof(char));
	outfilename = (char*)calloc(getfnsize(outputfile) + 50, sizeof(char));
	getfilename(inputfile, infilename);
	strcat(outfilename,"out_");
	strcat(outfilename, outputfile.filename);
	strcpy(outputfile.filename, outfilename);

	/* WAVEファイルからステレオの音データを入力する */
	stereo_wave_read(&org_pcm, infilename);
	printf("ファイルを読み込みました\n");
	/*for (n = 0; n < org_pcm.length; n++) {
		printf("org_pcm.sL[%d] = %lf\n", n, org_pcm.sL[n]);
	}*/
	printf("ファイル名：%s\n", infilename);
	printf("総データ数：%d[個]\n", (org_pcm.length) * 2);


	//ノイズ設定
	noise_index1 = getstereoindex(noise_time1, org_pcm.fs);
	noise_index2 = getstereoindex(noise_time2, org_pcm.fs);
	noise_size = noise_index2 - noise_index1;
	printf("\n----ノイズに関するデータ----\n");
	printf("ノイズ開始時刻：%lf, ノイズ終了時刻：%lf\n", noise_time1, noise_time2);
	printf("ノイズ開始データ = pcm[%d], ノイズ終了データ = pcm[%d]\n", noise_index1, noise_index2);
	setspcm(&noise_pcm, &org_pcm, noise_size / 2);
	getpcm(&org_pcm, &noise_pcm, noise_index1, noise_index2);
	/*for (n = 0; n < noise_size/2; n++) {
		printf("noise_pcm.sL[%d] = %lf\n", n, noise_pcm.sL[n]);
	}*/

	n_real = (double*)calloc(noise_size, sizeof(double));
	n_imag = (double*)calloc(noise_size, sizeof(double));

	//推定ノイズの初期化
	for (n = 0; n < noise_size; n++)
	{
		n_real[n] = 0.0;
		n_imag[n] = 0.0;
	}
	setnoise(n_real, n_imag,noise_pcm);
	/*for (n = 0; n < noise_size;n++) {
		printf("n_real[%d] = %lf\n",n,n_real[n]);
	}*/
	noiseSD_t = getSD(n_real, 0.0, noise_size);
	threshold = noiseSD_t * alpha;
	printf("平均：%lf\n", 0.0);
	printf("標準偏差：%lf\n", noiseSD_t);
	printf("閾値：%lf\n", threshold);

	spcmcpy(&out_pcm, &org_pcm);

	/* フレームの数 */
	number_of_frame = (org_pcm.length * 2 - L / 2) / (L / 2);
	printf("フレーム数：%d\n", number_of_frame);

	printf("\n----カウントに関するデータ----\n");
	sleepframe = getframe(sleeptime, L, org_pcm.fs);
	count = 0;
	printf("閾値：%lf\nカウント休止時間：%lf[s]\nカウント休止フレーム数：%d\n", threshold, sleeptime, sleepframe);


	/* メモリの確保 */
	x_real = (double*)calloc(N, sizeof(double));
	x_imag = (double*)calloc(N, sizeof(double));

	y_real = (double*)calloc(N, sizeof(double));
	y_imag = (double*)calloc(N, sizeof(double));


	printf("\n----解析開始----\n");
	//フレーム単位の周波数領域での処理
	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = (L / 2) * frame;

		/*データを0で初期化・入力*/
		for (n = 0; n < N; n++){
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
			y_real[n] = 0.0;
			y_imag[n] = 0.0;
		}
		for (n = 0; n < L / 2; n++){
			x_real[2 * n] = org_pcm.sL[(offset / 2) + n];
			x_real[(2 * n) + 1] = org_pcm.sR[(offset / 2) + n];
		}
		
		//データのセーブ
		for (n = 0; n < N; n++){
			y_real[n] = x_real[n];
			y_imag[n] = x_imag[n];
		}

	
		//音のカウント
		if (bflag == 1) {
			secondindex = getMaxindex(y_real, N);
			dy = representive_value - y_real[secondindex];
			if (dy > judgedy) {
				count++;
				temptime = getsecond(frame, L, org_pcm.fs);
				if (sleeptime > (temptime - counttime)) {
					sleeptime -= (temptime - counttime) * s;
				}
				else if (sleeptime < (temptime - counttime)) {
					sleeptime += (temptime - counttime) * s;
				}
				counttime = temptime;
				sleepframe = getframe(sleeptime, L, org_pcm.fs);
				printf("frame = %d, time = %lf, count = %d\n",frame-1, counttime, count);
				printf("カウント休止時間：%lf[s]\nカウント休止フレーム数：%d\n", sleeptime, sleepframe);
				cflag = sleepframe;

			}
			bflag = 0;
		}
		if (cflag == 0) {
			max = -1.0;
			min = 1.0;
			for (int i = 0; i < N; i++) {
				data = y_real[i];
				if (max < data) {
					max = data;
					maxindex = i;
				}
				if (min > data) {
					min = data;
					minindex = i;
				}
			}
			if (max > -min) {
				representive_value = max;
				representive_index = maxindex;
			}else {
				representive_value = min;
				representive_index = minindex;
			}
			if (representive_value > threshold) {
				printf("over:%lf[s]\n",getsecond(frame,L,org_pcm.fs));
				bflag = 1;
			}
		}else{
			cflag--;
		}

		/* フレームの連結 */
		for (n = 0; n < N / 2; n++)
		{
			if ((offset / 2) + n < org_pcm.length)
			{
				out_pcm.sL[(offset / 2) + n] = y_real[2 * n];
				out_pcm.sR[(offset / 2) + n] = y_real[(2 * n) + 1];
			}
		}
	}

	printf("解析時間： %lf\n", getsecond(number_of_frame, L, org_pcm.fs));
	printf("----解析終了----\n");


	/* WAVEファイルにモノラルの音データを出力する */
	getfilename(outputfile, outfilename);
	stereo_wave_write(&out_pcm, outfilename);
	printf("ファイルを出力しました\n");
	printf("ファイル名：%s\n", outfilename);

	/* メモリの解放 */
	free(x_imag);
	free(x_real);
	free(n_imag);
	free(n_real);
	free(outfilename);
	free(out_pcm.sL);
	free(out_pcm.sR);
	free(y_real);
	free(y_imag);
	free(infilename);
	free(org_pcm.sR);
	free(org_pcm.sL);

	printf("finish\n");

	return 0;
}
