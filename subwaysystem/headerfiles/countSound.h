#pragma once
#include "wave.h"
#include "dfr_filter.h"
#include "statistics.h"
#include "frame-time.h"

int countsound(double sound[], int soundlength, int framelength, double threshold) {
	

	return 0;
}

int judgewake() {
	int judge = 1;
	return 1;
}

//�P�t���[�����̉��f�[�^�A�t���[�����A������臒l��ݒ�
//�^����ꂽ�t���[�����̍ő剹��臒l�𒴂��邩����
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


//�U���X�y�N�g���Afft�|�C���g�A��͎��g���͈́A�������A臒l�����
int judgeFrequencyBySd(double *A, int fftpt, double fe1, double fe2, int divisionFe, int SamplesPerSec, double upperThreshold,double lowerThreshold) {

	int judge,section,sectionLength, frequencyOffset,divisionNumber,
	position,fe1_position,fe2_position,positionDistance;
	double *max,data,myu,sd;
	
	judge = 0;

	//�w�肵�����g�������̒l���܂܂��z��̃C���f�b�N�X
	fe1_position =(int) floor((fe1 * fftpt) / SamplesPerSec);		
	fe2_position =(int) floor((fe2 * fftpt) / SamplesPerSec);
	sectionLength =(int) floor((divisionFe * fftpt) / SamplesPerSec);

	//�s�b�^��sectionlength�Ŋ���؂��悤�ɒ���
	positionDistance = (fe2_position - fe1_position);
	divisionNumber = positionDistance / sectionLength;
	fe2_position -= positionDistance % sectionLength;
	max = (double*)calloc(divisionNumber, sizeof(double));

	for (section = 0; section < divisionNumber;section++) {
		frequencyOffset = section * sectionLength;
		max[section] = 0.0;
		for (position = 0; position < sectionLength; position++) {
				data = A[fe1_position + frequencyOffset + position];
				if (max[section] < data) {
					max[section] = data;
				}
			}
		}

	myu = getmyu(max, divisionNumber);
	sd = getSD(max,myu,divisionNumber);

	//if (sd < upperThreshold && sd > lowerThreshold) {
	//	printf("sd = %lf	", sd);
	//	judge = 1;
	//}
	if (myu < upperThreshold && myu > lowerThreshold) {
		printf("myu = %lf	", myu);
		judge = 1;
	}


	return judge;
}

int judgeFrequencyByAverageDif(double *A,double *average,int generation, int fftpt, double fe1, double fe2, int divisionFe, int SamplesPerSec, double magnification) {

	int judge, section, sectionLength, frequencyOffset, divisionNumber,
		position, fe1_position, fe2_position, positionDistance,
		g;
	double *max, data, myu,threshold;

	judge = 0;

	//�w�肵�����g�������̒l���܂܂��z��̃C���f�b�N�X
	fe1_position = (int)floor((fe1 * fftpt) / SamplesPerSec);
	fe2_position = (int)floor((fe2 * fftpt) / SamplesPerSec);
	sectionLength = (int)floor((divisionFe * fftpt) / SamplesPerSec);

	//�s�b�^��sectionlength�Ŋ���؂��悤�ɒ���
	positionDistance = (fe2_position - fe1_position);
	divisionNumber = positionDistance / sectionLength;
	fe2_position -= positionDistance % sectionLength;
	max = (double*)calloc(divisionNumber, sizeof(double));


	for (section = 0; section < divisionNumber; section++) {
		frequencyOffset = section * sectionLength;
		max[section] = 0.0;
		for (position = 0; position < sectionLength; position++) {
			data = A[fe1_position + frequencyOffset + position];
			if (max[section] < data) {
				max[section] = data;
			}
		}
	}

	myu = getmyu(max, divisionNumber);
	threshold = myu * magnification;

	g = 0;
	while (g < generation)
	{
		if (myu > average[g]*magnification) {
			printf("av:%f",average[g]);
			break;
		}
		g++;
	}
	if (g < generation) {
		printf("myu = %lf	", myu);
		judge = 1;
	}

	//average �̍X�V
	g = 0;
	while (g < generation - 1) {
	average[g] = average[g+1];
	g++;
	} 
	average[g] = myu;
	return judge;
}

int judgeFrequency(double *A, double *average, int generation, int fftpt, double fe1, double fe2, int divisionFe, int SamplesPerSec, double magnification, double upperThreshold, double lowerThreshold) {

	int judge, section, sectionLength, frequencyOffset, divisionNumber,
		position, fe1_position, fe2_position, positionDistance,
		g,flag;
	double *max, data, myu, threshold;

	judge = 0;

	//�w�肵�����g�������̒l���܂܂��z��̃C���f�b�N�X
	fe1_position = (int)floor((fe1 * fftpt) / SamplesPerSec);
	fe2_position = (int)floor((fe2 * fftpt) / SamplesPerSec);
	sectionLength = (int)floor((divisionFe * fftpt) / SamplesPerSec);

	//�s�b�^��sectionlength�Ŋ���؂��悤�ɒ���
	positionDistance = (fe2_position - fe1_position);
	divisionNumber = positionDistance / sectionLength;
	fe2_position -= positionDistance % sectionLength;
	max = (double*)calloc(divisionNumber, sizeof(double));


	for (section = 0; section < divisionNumber; section++) {
		frequencyOffset = section * sectionLength;
		max[section] = 0.0;
		for (position = 0; position < sectionLength; position++) {
			data = A[fe1_position + frequencyOffset + position];
			if (max[section] < data) {
				max[section] = data;
			}
		}
	}

	myu = getmyu(max, divisionNumber);

	g = 0;
	if (myu < upperThreshold && myu > lowerThreshold) {
		while (g < generation)
		{
			g++;
		}
		printf("myu = %lf	", myu);
		judge = 1;
	}

	//average �̍X�V
	g = 0;
	while (g < generation - 1) {
		average[g] = average[g + 1];
		g++;
	}
	average[g] = myu;
	return judge;
}
//������
double updateSleepTime(int wakeFlag ,double sleeptime,double s,double t,double before_counttime, int offset,int fs,int channel) {
	double temptime,newsleeptime;
	if (wakeFlag == 1 ) {
		temptime = getSecond(offset, fs,channel);
		if (sleeptime > (temptime - before_counttime)*t) {
			newsleeptime = sleeptime -(temptime - before_counttime) * s;
		}
		else if (sleeptime < (temptime - before_counttime)*t) {
			newsleeptime = sleeptime + (temptime - before_counttime) * s;
		}
		
	}
	else {
		newsleeptime = sleeptime;
	}
	return newsleeptime;
	
}

double updateValue(double value, double factor, double effetrate) {
	return factor * effetrate + value * (1 - effetrate);
}
