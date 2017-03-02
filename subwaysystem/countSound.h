#pragma once
#include "wave.h"
#include "dfr_filter.h"
#include "frame-time.h"
#include "statistics.h"
#define M 2000

int countsound(double sound[], int soundlength, int framelength, double threshold) {


	return 0;
}

int judgeSounnd(double sound[], int N, double threshold) {

	int count = 0;
	return count;
}

int CountTimes(int N, double sleeptime, double s, double t, int number_fe1, int number_fe2, double threshold, double threshold_2, double threshold_3, STEREO_PCM* org_pcm)
{
	double temptime, counttime, sum_all, data_all[M], *A, *A_all, temp, *w_real, *w_imag, *x_imag, *x_real, *y_imag, *y_real, *w, *data_ex, myu, sd, data_sd[M];
	int maxindex, minindex, representive_index, secondindex, bflag, cflag,sleepframe, count, frame, data, number_of_frame, i, j, k, m, n, offset, judge;
	bflag = 0;	//前フレームで閾値を超えたか判定
	cflag = 0;	//カウントをチェックするかどうか判定
	counttime = 0.0;	//前のカウント時刻
	count = 0;
	w = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	Hanning_window(w, N); /* ハニング窓 */

	x_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	w_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	w_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	A = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	A_all = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	number_of_frame = org_pcm->length / N;
	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = (N / 2) * frame;

		int maxindex, minindex, index;
		double max, min, data_ex, delta;

		for (n = 0; n < N; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
		}
		for (n = 0; n < N / 2; n++)
		{
			x_real[2 * n] = org_pcm->sL[offset + n];
			x_real[2 * n + 1] = org_pcm->sR[offset + n];
		}
		
		max = -1.0;
       	min = 1.0;
		

			for (i = 0; i < N; i++) 
			{
				data_ex = x_real[i];
				if (max < data_ex) {
					max = data_ex;
					maxindex = i;
				}
				if (min > data_ex) {
					min = data_ex;
					minindex = i;
				}
				
			}
			

		/*delta = max - min;*/
		if (max > -min) {
			delta = max;
			index = maxindex;
		}
		else {
			delta = -min;
			index = minindex;
		}

		if (delta > threshold_2) {
			judge = 1;
		}
		else{
			judge = 0;
		}
		/*printf("delta = %f\n", delta);*/
		

		/* w(m)のFFT */
		for (m = 0; m < N; m++)
		{
			w_real[m] = 0.0;
			w_imag[m] = 0.0;
		}
		for (m = 0; m < N; m++)
		{
			w_real[m] = w[m];
		}
		/*FFT(b_real, b_imag, N);*/

		/* 掛け合わせ */
		for (k = 0; k < N; k++)
		{
			y_real[k] = x_real[k] * w_real[k];
			y_imag[k] = x_imag[k] * w_imag[k];
		}
		FFT(y_real, y_imag, N);


		double A_all[M];
		for (k = 0; k < N; k++)
		{
			A_all[k] = sqrt(y_real[k] * y_real[k] + y_imag[k] * y_imag[k]);
		}

		/* 上位1/4の合計を利用した判定 */
		for (i = number_fe1; i < number_fe2 + 1; i++)
		{
			data_all[i] = A_all[i];
			for (j = i + 1; j < number_fe2 + 1; j++)
			{
				if (data_all[i] < data_all[j])
				{
					temp = data_all[i];
					data_all[i] = data_all[j];
					data_all[j] = temp;
				}
			}
		}

		sum_all = 0; /* 初期化 */
			for (i = number_fe1; i < (number_fe2 - number_fe1) / 4 + number_fe1 + 1; i++)
			{
				sum_all += data_all[i];
			}

			for (i = 0; i < N; i++)
			{
				data_sd[i] = A_all[i];
			}
			myu = getsoundmyu(data_sd, N);
			sd = getSD(data_sd, myu, N);

		if (cflag == 0) {
			
			if ((sum_all > threshold) && (judge == 1) && (sd<threshold_3)) {
				//printf("over:%lf[s]\n",getsecond(frame,L,org_pcm.fs));
				count++;
				temptime = getsecond(frame, N, org_pcm->fs);
				if (sleeptime > (temptime - counttime)*t) {
					sleeptime -= (temptime - counttime) * s;
				}
				else if (sleeptime < (temptime - counttime)*t) {
					sleeptime += (temptime - counttime) * s;
				}
				counttime = temptime;
				sleepframe = getframe(sleeptime, N, org_pcm->fs);
				/*printf("frame = %d, time = %lf, count = %d\n", frame - 1, counttime, count);
				printf("カウント休止時間：%lf[s]\nカウント休止フレーム数：%d\n", sleeptime, sleepframe);*/
				printf("time:%f    count = %d\n    delta = %f\n", (double)N * frame  / org_pcm->fs, count, delta);
				cflag = sleepframe;
				//bflag = 1;
			}
		}
		else {
			cflag--;
		}	
	}
	printf("回数は%d回です\n", count);
	return count;
}

int CountTimes_mono(int N, double sleeptime, double s, double t, int number_fe1, int number_fe2, double threshold, double threshold_2, double threshold_3, MONO_PCM* org_pcm)
{
	double temptime, counttime, sum_all, data_all[M], *A, *A_all, temp, *w_real, *w_imag, *x_imag, *x_real, *y_imag, *y_real, *w, *data_ex, data_sd[M], myu, sd;
	int maxindex, minindex, representive_index, secondindex, bflag, cflag, sleepframe, count, frame, data, number_of_frame, i, j, k, m, n, offset, judge;
	bflag = 0;	//前フレームで閾値を超えたか判定
	cflag = 0;	//カウントをチェックするかどうか判定
	counttime = 0.0;	//前のカウント時刻
	count = 0;
	w = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	Hanning_window(w, N); /* ハニング窓 */

	x_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	x_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	y_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	w_real = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	w_imag = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	A = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	A_all = (double*)calloc(N, sizeof(double)); /* メモリの確保 */
	number_of_frame = org_pcm->length / N;
	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = N * frame;

		int maxindex, minindex, index;
		double max, min, data_ex, delta;

		for (n = 0; n < N; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
		}
		for (n = 0; n < N; n++)
		{
			x_real[n] = org_pcm->s[offset + n];
		}

		max = -1.0;
		min = 1.0;


		for (i = 0; i < N; i++)
		{
			data_ex = x_real[i];
			if (max < data_ex) {
				max = data_ex;
				maxindex = i;
			}
			if (min > data_ex) {
				min = data_ex;
				minindex = i;
			}

		}


		/*delta = max - min;*/
		if (max > -min) {
			delta = max;
			index = maxindex;
		}
		else {
			delta = -min;
			index = minindex;
		}

		if (delta > threshold_2) {
			judge = 1;
		}
		else {
			judge = 0;
		}
		/*printf("delta = %f\n", delta);*/


		/* w(m)のFFT */
		for (m = 0; m < N; m++)
		{
			w_real[m] = 0.0;
			w_imag[m] = 0.0;
		}
		for (m = 0; m < N; m++)
		{
			w_real[m] = w[m];
		}
		/*FFT(b_real, b_imag, N);*/

		/* 掛け合わせ */
		for (k = 0; k < N; k++)
		{
			y_real[k] = x_real[k] * w_real[k];
			y_imag[k] = x_imag[k] * w_imag[k];
		}
		FFT(y_real, y_imag, N);


		double A_all[M];
		for (k = 0; k < N; k++)
		{
			A_all[k] = sqrt(y_real[k] * y_real[k] + y_imag[k] * y_imag[k]);
		}

		/* 上位1/4の合計を利用した判定 */
		for (i = number_fe1; i < number_fe2 + 1; i++)
		{
			data_all[i] = A_all[i];
			for (j = i + 1; j < number_fe2 + 1; j++)
			{
				if (data_all[i] < data_all[j])
				{
					temp = data_all[i];
					data_all[i] = data_all[j];
					data_all[j] = temp;
				}
			}
		}

		sum_all = 0; /* 初期化 */
		for (i = number_fe1; i < (number_fe2 - number_fe1) / 4 + number_fe1 + 1; i++)
		{
			sum_all += data_all[i];
		}

		for (i = 0; i < N; i++)
		{
			data_sd[i] = A_all[i];
		}
		myu = getsoundmyu(data_sd, N);
		sd = getSD(data_sd, myu, N);

		if (cflag == 0) {

			if ((sum_all>threshold) && (judge == 1) && (sd < threshold_3)) {
				//printf("over:%lf[s]\n",getsecond(frame,L,org_pcm.fs));
				count++;
				temptime = getsecond(frame, N, org_pcm->fs);
				if (sleeptime > (temptime - counttime)*t) {
					sleeptime -= (temptime - counttime) * s;
				}
				else if (sleeptime < (temptime - counttime)*t) {
					sleeptime += (temptime - counttime) * s;
				}
				counttime = temptime;
				sleepframe = getframe(sleeptime, N, org_pcm->fs);
				/*printf("frame = %d, time = %lf, count = %d\n", frame - 1, counttime, count);
				printf("カウント休止時間：%lf[s]\nカウント休止フレーム数：%d\n", sleeptime, sleepframe);*/
				printf("time:%f    count = %d\n    標準偏差: %f\n", (double)N * frame / org_pcm->fs, count, sd);
				cflag = sleepframe;
				//bflag = 1;
			}
		}
		else {
			cflag--;
		}
	}
	printf("回数は%d回です\n", count);
	return count;
}