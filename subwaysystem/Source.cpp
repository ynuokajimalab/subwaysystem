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
	STEREO_PCM org_pcm, out_pcm;
	FILENAME inputfile,outputfile;
	int n, N, L, offset, frame, number_of_frame;
	int cflag, count;
	double *x_real, *x_imag, *y_real,*y_imag;
	double threshold, sleeptime, alpha, counttime, s,t;
	char *infilename, *outfilename;


	//入力ファイルのデータ
	char orgfile[] = "short2";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//フレームの長さ
	N = 1024;
	L = 512;  //実際のデータサイズ
	threshold = 0;
	sleeptime = 0.5;
	alpha = 4.4;
	s = 0.06;
	t = 2.0/3.0;


	//データの初期化
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
	printf("ファイル名：%s\n", infilename);
	printf("総データ数：%d[個]\n", (org_pcm.length) * 2);

	spcmcpy(&out_pcm, &org_pcm);

	/* フレームの数 */
	number_of_frame = (org_pcm.length * 2 - L / 2) / (L / 2);
	printf("フレーム数：%d\n", number_of_frame);

	count = 0;
	printf("閾値：%lf\nカウント休止時間：%lf[s]\n", threshold, sleeptime);


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

		if (cflag == 0) {
			if (judgeSoundPower(x_real,N,threshold) == 1) {
				if (true) {
				count++;
				//updateSleepTime(sleeptime, counttime, offset, org_pcm.fs);
				}
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
	printf("総カウント数： %d\n",count);
	printf("----解析終了----\n");


	/* WAVEファイルにモノラルの音データを出力する */
	getfilename(outputfile, outfilename);
	stereo_wave_write(&out_pcm, outfilename);
	printf("ファイルを出力しました\n");
	printf("ファイル名：%s\n", outfilename);

	/* メモリの解放 */
	free(x_imag);
	free(x_real);
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
