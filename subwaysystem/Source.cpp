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


	//���̓t�@�C���̃f�[�^
	char orgfile[] = "short2";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//�t���[���̒���
	N = 1024;
	L = 512;  //���ۂ̃f�[�^�T�C�Y
	threshold = 0;
	sleeptime = 0.5;
	alpha = 4.4;
	s = 0.06;
	t = 2.0/3.0;


	//�f�[�^�̏�����
	cflag = 0;	//�J�E���g���`�F�b�N���邩�ǂ�������
	counttime = 0.0;	//�O�̃J�E���g����

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
	number_of_frame = (org_pcm.length * 2 - L / 2) / (L / 2);
	printf("�t���[�����F%d\n", number_of_frame);

	count = 0;
	printf("臒l�F%lf\n�J�E���g�x�~���ԁF%lf[s]\n", threshold, sleeptime);


	/* �������̊m�� */
	x_real = (double*)calloc(N, sizeof(double));
	x_imag = (double*)calloc(N, sizeof(double));

	y_real = (double*)calloc(N, sizeof(double));
	y_imag = (double*)calloc(N, sizeof(double));


	printf("\n----��͊J�n----\n");
	//�t���[���P�ʂ̎��g���̈�ł̏���
	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = (L / 2) * frame;

		/*�f�[�^��0�ŏ������E����*/
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
