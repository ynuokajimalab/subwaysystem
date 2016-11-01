#pragma warning(disable:4996)

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "./headerfiles/wave.h"
#include "./headerfiles/window_function.h"
#include "./headerfiles/sinc.h"
#include "./headerfiles/fir_filter.h"
#include "./headerfiles/fft.h"
#include "./headerfiles/countSound.h"
#include "./headerfiles/make_datfile.h"
#include "./headerfiles/dfr_filter.h"


int main(void)
{
	STEREO_PCM pcm0,pcm0_edge, pcm1,pcm1_edge;
	int n, m, k, J, L, N,ife1,ife2, offset, frame, number_of_frame;
	double fe1, fe2, delta,xedge_max,yedge_max, tempLmax, tempRmax, *b, *w, *b_real, *b_imag, *x_real, *x_imag, *x_edge, *y_real, *y_imag, *y_edge;
	char *inputfile, *outputfile,*outputedgefile, *frequency;
	size_t infsize, outfsize;

	//���̓t�@�C���̐ݒ�
	char filename[] = "short1";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	infsize = strlen(directory)+ strlen(filename)+ strlen(filetype)+1;
	inputfile =(char*)calloc(infsize ,sizeof(char));
	strcpy(inputfile, directory);
	strcat(inputfile ,filename);
	strcat(inputfile, filetype);


	stereo_wave_read(&pcm0,inputfile); /* WAVE�t�@�C������X�e���I�̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	pcm0_edge.fs = pcm0.fs; /* �W�{�����g�� */
	pcm0_edge.bits = pcm0.bits; /* �ʎq�����x */
	pcm0_edge.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm0_edge.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm0_edge.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	pcm1_edge.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1_edge.bits = pcm0.bits; /* �ʎq�����x */
	pcm1_edge.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1_edge.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm1_edge.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	fe1 = 1200.0 / (2 * pcm0.fs); /* �G�b�W���g��1 */
	fe2 = 2000.0 / (2 * pcm0.fs); /* �G�b�W���g��2 */
	delta = 800.0 / (2 * pcm0.fs); /* �J�ڑш敝 */
	ife1 = (int)(fe1 * 2 * pcm0.fs);
	ife2 = (int)(fe2 * 2 * pcm0.fs);
	//5-5��
	frequency = (char*)calloc(13, sizeof(char));
	sprintf(frequency,"_%d-%d",ife1,ife2);

	J = (int)(3.1 / delta + 0.5) - 1; /* �x����̐� */
	if (J % 2 == 1)
	{
		J++; /* J+1����ɂȂ�悤�ɒ������� */
	}

	b = (double*) calloc((J + 1), sizeof(double)); /* �������̊m�� */
	w = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */

	Hanning_window(w, (J + 1)); /* �n�j���O�� */

	FIR_BPF(fe1,fe2, J, b, w); /* FIR�t�B���^�̐݌v */

	L = 256; /* �t���[���̒��� */
	N = 512; /* DFT�̃T�C�Y */

	yedge_max = 0; //edge�f�[�^�̍ő�l

	number_of_frame = pcm0.length*2 / L; /* �t���[���̐� */

	 /* �������̊m�� */
	b_real = (double*)calloc(N, sizeof(double)); 
	b_imag = (double*)calloc(N, sizeof(double)); 
	x_real = (double*)calloc(N, sizeof(double));
	x_imag = (double*)calloc(N, sizeof(double));
	x_edge = (double*)calloc(N, sizeof(double));
	y_real = (double*)calloc(N, sizeof(double));
	y_imag = (double*)calloc(N, sizeof(double));
	y_edge = (double*)calloc(N, sizeof(double));

	/*�G�b�W��������*/
	laplacian(pcm0.sL, pcm0_edge.sL, pcm0_edge.length);
	laplacian(pcm0.sR, pcm0_edge.sR, pcm0_edge.length);
	//edgestress(y_real, y_edge, N,enhancerate);
	tempLmax = getMax(pcm0_edge.sL, pcm0_edge.length);
	tempRmax = getMax(pcm0_edge.sR, pcm0_edge.length);
	if (tempLmax < tempRmax) {
		xedge_max = tempRmax;
	}
	else {
		xedge_max = tempLmax;
	}
	//memcpy(pcm0.sR, pcm0_edge.sR, pcm0_edge.length);
	//memcpy(pcm0.sR, pcm0_edge.sR, pcm0_edge.length);

	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;

		/* x(n)��FFT */
		for (n = 0; n < N; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
			x_edge[n] = 0.0;
		}
		for (n = 0; n < L/2; n++)
		{
			x_real[2 * n] = pcm0.sL[(offset/2) + n];
			x_real[(2 * n) + 1] = pcm0.sR[(offset / 2) + n];
		}

		/*�G�b�W��������*/
		//laplacian(x_real, x_edge, N);
		//edgestress(y_real, y_edge, N,enhancerate);

	    //memcpy(x_real,x_edge,N);

		FFT(x_real, x_imag, N);

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
			y_real[k] = x_real[k] * b_real[k] - x_imag[k] * b_imag[k];
			y_imag[k] = x_imag[k] * b_real[k] + x_real[k] * b_imag[k];

			y_edge[k] = 0;
		}

		IFFT(y_real, y_imag, N);

		/*�G�b�W��������*/
		laplacian(y_real, y_edge, N);
		//edgestress(y_real, y_edge, N,enhancerate);


		/* �t�B���^�����O���ʂ̘A�� */
		for (n = 0; n < L * 2/2; n++)
		{
			if ((offset/2) + n < pcm1.length)
			{	
				pcm1.sL[(offset / 2) + n] += y_real[2 * n];
				pcm1.sR[(offset / 2) + n] += y_real[(2 * n) + 1];

				pcm1_edge.sL[(offset/2) + n] += y_edge[2 * n];			
				pcm1_edge.sR[(offset/2) + n] += y_edge[(2*n)+1];
			}
		}
	}

	/*�G�b�W��������*/
	laplacian(pcm1.sL, pcm1_edge.sL, pcm1_edge.length);
	laplacian(pcm1.sR, pcm1_edge.sR, pcm1_edge.length);
	//edgestress(y_real, y_edge, N,enhancerate);


	tempLmax = getMax(pcm1_edge.sL,pcm1_edge.length);
	tempRmax = getMax(pcm1_edge.sR, pcm1_edge.length);

	if (tempLmax < tempRmax) {
		yedge_max = tempRmax;
	}else{
		yedge_max = tempLmax;
	}

	
	//�������x���̐��K���A�X�y�N�g���͂Ȃ܂��Č�����
	normalize_array(pcm1_edge.sL, pcm1_edge.length, yedge_max);
	normalize_array(pcm1_edge.sR, pcm1_edge.length, yedge_max);

	//�t�@�C����
	outfsize = strlen(directory) + 4 + strlen(filename) + strlen(frequency) +strlen(filetype) + 1;
	outputfile = (char*)calloc(outfsize,sizeof(char));
	outputedgefile = (char*)calloc(outfsize+5, sizeof(char));
	strcpy(outputfile, directory);
	strcat(outputfile, "out_");
	strcpy(outputedgefile, outputfile);
	strcat(outputedgefile,"edge_");
	strcat(outputfile, filename);
	strcat(outputfile, frequency);
	strcat(outputfile,filetype);
	strcat(outputedgefile, filename);
	strcat(outputedgefile, frequency);
	strcat(outputedgefile, filetype);

	stereo_wave_write(&pcm1, outputfile); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
	stereo_wave_write(&pcm1_edge, outputedgefile); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */

	 /* �������̉�� */
	free(frequency);
	free(inputfile);
	free(outputfile);
	free(outputedgefile);
	free(pcm0.sL);
	free(pcm0.sR);
	free(pcm1.sL); 
	free(pcm1.sR);
	free(pcm0_edge.sL);
	free(pcm0_edge.sR);
	free(pcm1_edge.sL); 
	free(pcm1_edge.sR); 
	free(b); 
	free(w);
	free(b_real);
	free(b_imag);
	free(x_real);
	free(x_imag);
	free(x_edge);
	free(y_real);
	free(y_imag);
	free(y_edge);


	//printf("counts =%d\n count = %d\n",counts,count);
	//printf("max = %lf\n ",max);
	printf("finish");

	return 0;
}
