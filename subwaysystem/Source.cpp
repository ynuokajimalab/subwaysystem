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


	//���̓t�@�C���̃f�[�^
	char orgfile[] = "short3";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//�t���[���̒���
	N = 1024;
	L = 512;  //���ۂ̃f�[�^�T�C�Y
	//�m�C�Y�̃f�[�^
	noise_time1 = 0.0;
	noise_time2 = 7.0;
	sleeptime = 0.3;
	judgedy = 0.04;
	alpha = 4.9;
	s = 0.02;


	//�f�[�^�̏�����
	bflag = 0;	//�O�t���[����臒l�𒴂���������
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
	/*for (n = 0; n < org_pcm.length; n++) {
		printf("org_pcm.sL[%d] = %lf\n", n, org_pcm.sL[n]);
	}*/
	printf("�t�@�C�����F%s\n", infilename);
	printf("���f�[�^���F%d[��]\n", (org_pcm.length) * 2);


	//�m�C�Y�ݒ�
	noise_index1 = getstereoindex(noise_time1, org_pcm.fs);
	noise_index2 = getstereoindex(noise_time2, org_pcm.fs);
	noise_size = noise_index2 - noise_index1;
	printf("\n----�m�C�Y�Ɋւ���f�[�^----\n");
	printf("�m�C�Y�J�n�����F%lf, �m�C�Y�I�������F%lf\n", noise_time1, noise_time2);
	printf("�m�C�Y�J�n�f�[�^ = pcm[%d], �m�C�Y�I���f�[�^ = pcm[%d]\n", noise_index1, noise_index2);
	setspcm(&noise_pcm, &org_pcm, noise_size / 2);
	getpcm(&org_pcm, &noise_pcm, noise_index1, noise_index2);
	/*for (n = 0; n < noise_size/2; n++) {
		printf("noise_pcm.sL[%d] = %lf\n", n, noise_pcm.sL[n]);
	}*/

	n_real = (double*)calloc(noise_size, sizeof(double));
	n_imag = (double*)calloc(noise_size, sizeof(double));

	//����m�C�Y�̏�����
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
	printf("���ρF%lf\n", 0.0);
	printf("�W���΍��F%lf\n", noiseSD_t);
	printf("臒l�F%lf\n", threshold);

	spcmcpy(&out_pcm, &org_pcm);

	/* �t���[���̐� */
	number_of_frame = (org_pcm.length * 2 - L / 2) / (L / 2);
	printf("�t���[�����F%d\n", number_of_frame);

	printf("\n----�J�E���g�Ɋւ���f�[�^----\n");
	sleepframe = getframe(sleeptime, L, org_pcm.fs);
	count = 0;
	printf("臒l�F%lf\n�J�E���g�x�~���ԁF%lf[s]\n�J�E���g�x�~�t���[�����F%d\n", threshold, sleeptime, sleepframe);


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
		
		//�f�[�^�̃Z�[�u
		for (n = 0; n < N; n++){
			y_real[n] = x_real[n];
			y_imag[n] = x_imag[n];
		}

	
		//���̃J�E���g
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
				printf("�J�E���g�x�~���ԁF%lf[s]\n�J�E���g�x�~�t���[�����F%d\n", sleeptime, sleepframe);
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

	printf("��͎��ԁF %lf\n", getsecond(number_of_frame, L, org_pcm.fs));
	printf("----��͏I��----\n");


	/* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
	getfilename(outputfile, outfilename);
	stereo_wave_write(&out_pcm, outfilename);
	printf("�t�@�C�����o�͂��܂���\n");
	printf("�t�@�C�����F%s\n", outfilename);

	/* �������̉�� */
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
