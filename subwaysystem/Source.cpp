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
#include "./headerfiles/fft.h"
#include "./headerfiles/countSound.h"
#include "./headerfiles/make_datfile.h"
#include "./headerfiles/dfr_filter.h"
#define OUT_FILE_NUM 4

int main(void)
{
	STEREO_PCM pcm0, bflap_pcm0, pcm1[OUT_FILE_NUM];
	IN_FILE inputfile;
	OUT_FILE outputfiles[OUT_FILE_NUM];
	int n, m, k, J, L, N,ife1,ife2,idelta, offset, frame, number_of_frame;
	double fe1, fe2, delta,xedge_max,yedge_max[OUT_FILE_NUM / 2], tempLmax, tempRmax, *b, *w, *b_real, *b_imag, *x_real, *x_imag, *x_lpre, *x_lpim, *y_real[OUT_FILE_NUM], *y_imag[OUT_FILE_NUM];
	char *infilename, *outfilenames[OUT_FILE_NUM];


	//���̓t�@�C���̃f�[�^
	char file[] = "thelatest";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	//BP�t�B���^�̃f�[�^
	ife1 = 2500;	/* �G�b�W���g��1 */
	ife2 = 3500;	/* �G�b�W���g��2 */
	idelta = 900;	/* �J�ڑш敝 */
	//DFT�̃f�[�^
	L = 256; /* �t���[���̒��� */
	N = 512; /* DFT�̃T�C�Y */


	//�t�@�C���̏���������
	init_infiles(&inputfile,directory,file,filetype);
	infilename =(char*)calloc(getinfsize(&inputfile) ,sizeof(char));
	getinfilename(&inputfile,infilename);
	init_outfiles(&inputfile,outputfiles,ife1,ife2,idelta);
	for (int i = 0; i < OUT_FILE_NUM; i++) {
		outfilenames[i] = (char*)calloc(getoutfsize(&outputfiles[i]), sizeof(char));
		getoutfilename(&outputfiles[i],outfilenames[i]);
	}

	stereo_wave_read(&pcm0,infilename); /* WAVE�t�@�C������X�e���I�̉��f�[�^����͂��� */

	spcmcpy(&bflap_pcm0,&pcm0);
	for (int i = 0; i < OUT_FILE_NUM;i++) {
		spcmcpy(&pcm1[i],&pcm0);
	}

	//BP�t�B���^�������
	fe1 = (double)ife1 / (2 * pcm0.fs); /* �G�b�W���g��1 */
	fe2 = (double)ife2 / (2 * pcm0.fs); /* �G�b�W���g��2 */
	delta = (double)idelta / (2 * pcm0.fs); /* �J�ڑш敝 */
	J = (int)(3.1 / delta + 0.5) - 1; /* �x����̐� */
	if (J % 2 == 1)
	{
		J++; /* J+1����ɂȂ�悤�ɒ������� */
	}
	b = (double*) calloc((J + 1), sizeof(double)); /* �������̊m�� */
	w = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */
	Hanning_window(w, (J + 1)); /* �n�j���O�� */
	FIR_BPF(fe1,fe2, J, b, w); /* FIR�t�B���^�̐݌v */


	yedge_max[0] = 0; //edge�f�[�^�̍ő�l
	yedge_max[1] = 0; //edge�f�[�^�̍ő�l

	number_of_frame = pcm0.length*2 / L; /* �t���[���̐� */

	 /* �������̊m�� */
	b_real = (double*)calloc(N, sizeof(double)); 
	b_imag = (double*)calloc(N, sizeof(double)); 
	x_real = (double*)calloc(N, sizeof(double));
	x_imag = (double*)calloc(N, sizeof(double));
	x_lpre = (double*)calloc(N, sizeof(double));
	x_lpim = (double*)calloc(N, sizeof(double));
	for (int i = 0; i < OUT_FILE_NUM; i++)
	{
		y_real[i] = (double*)calloc(N, sizeof(double));
		y_imag[i] = (double*)calloc(N, sizeof(double));
	}

	/*�G�b�W��������*/
	laplacian(pcm0.sL, bflap_pcm0.sL, bflap_pcm0.length);
	laplacian(pcm0.sR, bflap_pcm0.sR, bflap_pcm0.length);
	//tempLmax = getMax(bflap_pcm0.sL, bflap_pcm0.length);
	//tempRmax = getMax(bflap_pcm0.sR, bflap_pcm0.length);
	tempLmax = getbeforeMax(bflap_pcm0.sL, bflap_pcm0.length);
	tempRmax = getbeforeMax(bflap_pcm0.sR, bflap_pcm0.length);
	if (tempLmax < tempRmax) {
		xedge_max = tempRmax;
	}
	else {
		xedge_max = tempLmax;
	}

	//�������x���̐��K���A�X�y�N�g���͂Ȃ܂��Č�����
	normalize_array(bflap_pcm0.sL, bflap_pcm0.length, xedge_max);
	normalize_array(bflap_pcm0.sR, bflap_pcm0.length, xedge_max);


	//�t���[���P�ʂ̎��g���̈�ł̏���
	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;

		/* x(n)��FFT */
		for (n = 0; n < N; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
			x_lpre[n] = 0.0;
			x_lpim[n] = 0.0;
			for (int i = 0; i < OUT_FILE_NUM;i++) {
				y_real[i][n] = 0;
				y_imag[i][n] = 0;
			}
		}
		for (n = 0; n < L/2; n++)
		{
			x_real[2 * n] = pcm0.sL[(offset/2) + n];
			x_real[(2 * n) + 1] = pcm0.sR[(offset / 2) + n];
			x_lpre[2 * n] = bflap_pcm0.sL[(offset / 2) + n];
			x_lpre[(2 * n) + 1] = bflap_pcm0.sR[(offset / 2) + n];
		}

		/*�G�b�W��������*/
		//laplacian(x_real, x_edge, N);

		FFT(x_real, x_imag, N);
		FFT(x_lpre, x_lpim, N);

		/* b(m)��FFT */
		for (m = 0; m < N; m++)
		{
			b_real[m] = 0.0;
			b_imag[m] = 0.0;
		}
		for (m = 0; m <= J; m++)
		{
			b_real[m] = b[m];
		}
		FFT(b_real, b_imag, N);

		/* �t�B���^�����O */
		for (k = 0; k < N; k++)
		{
			y_real[0][k] = x_real[k] * b_real[k] - x_imag[k] * b_imag[k];
			y_imag[0][k] = x_imag[k] * b_real[k] + x_real[k] * b_imag[k];
			y_real[2][k] = x_lpre[k] * b_real[k] - x_lpim[k] * b_imag[k];
			y_imag[2][k] = x_lpim[k] * b_real[k] + x_lpre[k] * b_imag[k];
		}

		IFFT(y_real[0], y_imag[0], N);
		IFFT(y_real[2], y_imag[2], N);

		/*�G�b�W��������*/
		//laplacian(y_real, y_edge, N);

		/* �t�B���^�����O���ʂ̘A�� */
		for (n = 0; n < L * 2/2; n++)
		{
			if ((offset/2) + n < pcm0.length)
			{	
				pcm1[0].sL[(offset / 2) + n] += y_real[0][2 * n];
				pcm1[0].sR[(offset / 2) + n] += y_real[0][(2 * n) + 1];

				pcm1[2].sL[(offset/2) + n] += y_real[2][2 * n];			
				pcm1[2].sR[(offset/2) + n] += y_real[2][(2 * n)+1];
			}
		}
	}

	/*�G�b�W��������*/
	laplacian(pcm1[0].sL, pcm1[1].sL, pcm1[1].length);
	laplacian(pcm1[0].sR, pcm1[1].sR, pcm1[1].length);
	laplacian(pcm1[2].sL, pcm1[3].sL, pcm1[3].length);
	laplacian(pcm1[2].sR, pcm1[3].sR, pcm1[3].length);
	
	//�������x���̐��K���A�X�y�N�g���͂Ȃ܂��Č�����
	for (int i = 0; i < OUT_FILE_NUM / 2;i++) {
		tempLmax = getMax(pcm1[2*i +1].sL, pcm1[2*i +1].length);
		tempRmax = getMax(pcm1[2*i +1].sR, pcm1[2*i +1].length);
		if (tempLmax < tempRmax) {
			yedge_max[i] = tempRmax;
		}else{
			yedge_max[i] = tempLmax;
		}
	//normalize_array(pcm1[2 * i + 1].sL, pcm1[2 * i + 1].length, yedge_max[i]);
	//normalize_array(pcm1[2 * i + 1].sR, pcm1[2 * i + 1].length, yedge_max[i]);
	}

	for (int i = 0; i < OUT_FILE_NUM;i++) {
	stereo_wave_write(&pcm1[i], outfilenames[i]); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
	}
	

	 /* �������̉�� */
	free(infilename);
	free(pcm0.sL);
	free(pcm0.sR);
	free(bflap_pcm0.sL);
	free(bflap_pcm0.sR);
	free(b); 
	free(w);
	free(b_real);
	free(b_imag);
	free(x_real);
	free(x_imag);
	free(x_lpre);
	free(x_lpim);
	for (int i = 0; i < OUT_FILE_NUM; i++)
	{
	free(outfilenames[i]);
	free(pcm1[i].sL); 
	free(pcm1[i].sR);
	free(y_real[i]);
	free(y_imag[i]);
	}


	//printf("counts =%d\n count = %d\n",counts,count);
	//printf("max = %lf\n ",max);
	printf("finish");

	return 0;
}
