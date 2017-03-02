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
	int n,k,g, N, offset, frame, number_of_frame;
	int wakeflag, count,
		lowFrequency,highFrequency,divisionFrequency,generation;
	double *x_real, *x_imag, *y_real,*y_imag,*wN, *A,*average;
	double thresholdOfPower,upperThresholdOfFrequency, lowerThresholdOfFrequency, lowerThresholdRate,upperThresholdRate,
		sleeptime, alpha, counttime, s,t,
		frequencyDistance,magnification,effectRate;
	char *infilename, *outfilename;


	//���̓t�@�C���̃f�[�^
	char orgfile[] = "short3";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	int channel = 2;
	//�t���[���̒���
	N = 2048;
	//�������x������
	thresholdOfPower = 0;
	//���g������
	highFrequency = 3000;
	lowFrequency = 500;
	divisionFrequency = 50;
	lowerThresholdOfFrequency = 3.5;	//short1: myu > 3.5�t��
	upperThresholdOfFrequency = 10.0;
	lowerThresholdRate = 1.3;
	upperThresholdRate = 3.5;
	effectRate = 0.02;

	generation = 3;
	magnification = 0.0;

	sleeptime = 0.5;
	alpha = 4.4;
	s = 0.06;
	t = 2.0/3.0;


	//�f�[�^�̏�����
	wakeflag = 1;	//�J�E���g���`�F�b�N���邩�ǂ�������
	counttime = 0.0;	//�O�̃J�E���g����
	average  = (double*)calloc(generation, sizeof(double));
	for (g = 0; g < generation;g++) {
		average[g] = 100.0;
	}


	//�t�@�C���̏���������
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

	/* WAVE�t�@�C������X�e���I�̉��f�[�^����͂��� */
	stereo_wave_read(&org_pcm, infilename);
	printf("�t�@�C����ǂݍ��݂܂���\n");
	printf("�t�@�C�����F%s\n", infilename);
	printf("���f�[�^���F%d[��]\n", (org_pcm.length) * 2);

	spcmcpy(&out_pcm, &org_pcm);
	/* �t���[���̐� */
	number_of_frame = org_pcm.length*2/N;
	frequencyDistance = (double)org_pcm.fs/(double)(2*N);
	printf("�t���[�����F%d	�t���[�����F%f[s]\n", number_of_frame,getSecond(N,org_pcm.fs,2));
	printf("���g������\�F%f[Hz]\n", frequencyDistance);
	printf("�O�f�[�^�e���W���F%f	��臒l�������G%f	��臒l�������F%f\n",effectRate,lowerThresholdRate,upperThresholdRate);
	//printf("����F%d	臒l�ω����F%lf\n",generation,magnification);


	printf("---------���g������---------\n");
	printf("�͈́F%d�`%d[Hz]�A���%d[Hz]\n",lowFrequency,highFrequency,divisionFrequency);
	printf("臒l�Fsd[%f�`%f]�Ŕ���\n",lowerThresholdOfFrequency,upperThresholdOfFrequency);

	printf("---------�������x������---------\n");
	printf("臒l�F%f�Ŕ���\n", thresholdOfPower);

	//�T�C�YN�̑��֐��쐬
	wN = (double*)calloc(N, sizeof(double));
	Hanning_window(wN, N);

	count = 0;

	/* �������̊m�� */
	x_real = (double*)calloc(N, sizeof(double));
	x_imag = (double*)calloc(N, sizeof(double));
	y_real = (double*)calloc(N, sizeof(double));
	y_imag = (double*)calloc(N, sizeof(double));
	A = (double*)calloc(N, sizeof(double));

	printf("\n----��͊J�n----\n");
	//�t���[���P�ʂ̎��g���̈�ł̏���
	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = N * frame;

		/*�f�[�^��0�ŏ������E����*/
		for (n = 0; n < N; n++){
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
			y_real[n] = 0.0;
			y_imag[n] = 0.0;
		}
		for (n = 0; n < N / 2; n++){
			x_real[2 * n] = org_pcm.sL[(offset / 2) + n];
			x_real[(2 * n) + 1] = org_pcm.sR[(offset / 2) + n];
		}
		for(n = 0; n < N; n++) {
			y_real[n] = wN[n] * x_real[n];
		}

		if ((wakeflag = judgewake()) == 1) {
			if(true)//if (judgeSoundPower(x_real,N,thresholdOfPower) == 1) 
			{
				FFT(y_real, y_imag, N);
				for (k = 0; k < N; k++) {
					A[k] = sqrt(y_real[k] * y_real[k] + y_imag[k] * y_imag[k]);
				}
				if (judgeFrequency(A,average,generation,N,lowFrequency ,highFrequency,divisionFrequency, org_pcm.fs, magnification, upperThresholdOfFrequency,lowerThresholdOfFrequency)==1)
				{
					printf("time = %f	lowerThreshold = %f, upperThreshold = %f\n", getSecond(offset,org_pcm.fs,channel),lowerThresholdOfFrequency,upperThresholdOfFrequency);
					count++;
				}
			}
		}
		lowerThresholdOfFrequency = updateValue(lowerThresholdOfFrequency,lowerThresholdRate* average[generation - 1],effectRate);
		upperThresholdOfFrequency = updateValue(upperThresholdOfFrequency, upperThresholdRate* average[generation - 1], effectRate);

				//sleeptime = updateSleepTime(wakeflag,sleeptime, counttime,s,t, offset, org_pcm.fs);

		/* �t���[���̘A�� */
		for (n = 0; n < N / 2; n++)
		{
			if ((offset / 2) + n < org_pcm.length)
			{
				out_pcm.sL[(offset / 2) + n] = y_real[2 * n];
				out_pcm.sR[(offset / 2) + n] = y_real[(2 * n) + 1];
			}
		}
	}
	printf("���J�E���g���F %d\n",count);
	printf("----��͏I��----\n");


	/* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
	getfilename(outputfile, outfilename);
	stereo_wave_write(&out_pcm, outfilename);
	printf("�t�@�C�����o�͂��܂���\n");
	printf("�t�@�C�����F%s\n", outfilename);

	/* �������̉�� */
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
