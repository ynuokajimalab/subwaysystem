#pragma once
#include "wave.h"
#include "dfr_filter.h"
#include "statistics.h"
#include "frame-time.h"

int countsound(double sound[], int soundlength, int framelength, double threshold) {
	

	return 0;
}

//１フレーム分の音データ、フレーム長、音圧の閾値を設定
//与えられたフレーム内の最大音が閾値を超えるか判定
int judgeSoundPower(double *x,int fftpt,double thresholdOfPower) {

	int n,judge = 0;
	double max, min,data,maxpower;

	max = -1.0;
	min = 1.0;
	for (n = 0; n < fftpt; n++) {
		data = x[n];
		if (max < data) {
			max = data;
		}
		if (min > data) {
			min = data;
		}
	}
	if (max > -min) {
		maxpower= max;
	}
	else {
		maxpower = -min;
	}
	if (maxpower > thresholdOfPower) {
		judge = 1;
	}

	return judge;
}

//振幅スペクトル、fftポイント、解析周波数範囲、分割数、閾値を入力
int judgeFrequencyBySd(double *A, int fftpt, double fe1, double fe2, int divisionFe, int SamplesPerSec, double upperThreshold,double lowerThreshold) {

	int judge,section,sectionLength, frequencyOffset,divisionNumber,
	position,fe1_position,fe2_position,positionDistance;
	double *max,data,myu,sd;
	
	judge = 0;

	//指定した周波数成分の値が含まれる配列のインデックス
	fe1_position =(int) floor((fe1 * fftpt) / SamplesPerSec);		
	fe2_position =(int) floor((fe2 * fftpt) / SamplesPerSec);
	sectionLength =(int) floor((divisionFe * fftpt) / SamplesPerSec);

	//ピッタリsectionlengthで割り切れるように調節
	positionDistance = (fe2_position - fe1_position);
	divisionNumber = positionDistance / sectionLength;
	fe2_position -= positionDistance % sectionLength;
	max = (double*)calloc(divisionNumber, sizeof(double));

	printf("divisionNumber = %d, sectionLength = %d",divisionNumber,sectionLength);

	for (section = 0; section < divisionNumber;section++) {
		frequencyOffset = section * sectionLength;
		max[section] = 0.0;
		for (position = 0; position < sectionLength; position++) {
				data = A[frequencyOffset + position];
				if (max[section] < data) {
					max[section] = data;
				}
			}
		}

	myu = getmyu(max, divisionNumber);
	sd = getSD(max,myu,divisionNumber);

	if (sd < upperThreshold && sd > lowerThreshold) {
		judge = 1;
	}

	return judge;
}

double updateSleepTime(double sleeptime,double s,double t,double before_counttime, int offset,int fs) {
	double temptime;

	temptime = getSecond(offset,fs);
	if (sleeptime > (temptime - before_counttime)*t) {
		sleeptime -= (temptime - before_counttime) * s;
	}
	else if (sleeptime < (temptime - before_counttime)*t) {
		sleeptime += (temptime - before_counttime) * s;
	}
	return sleeptime;
}
