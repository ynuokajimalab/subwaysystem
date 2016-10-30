#pragma warning(disable:4996)

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "./hederfiles/wave.h"
#include "./hederfiles/window_function.h"
#include "./hederfiles/sinc.h"
#include "./hederfiles/fir_filter.h"
#include "./hederfiles/fft.h"
#include "./hederfiles/countSound.h"
#include "./hederfiles/make_datfile.h"
#include "./hederfiles/dfr_filter.h"


int main(void)
{
	STEREO_PCM pcm0,pcm0_edge, pcm1,pcm1_edge;
	int n, m, k, J, L, N, offset, frame, number_of_frame;
	double fe1, fe2, delta, temp,xedge_max,yedge_max, tempLmax, tempRmax, *b, *w, *b_real, *b_imag, *x_real, *x_imag, *y_real, *y_imag, *y_adge;

	stereo_wave_read(&pcm0,"./wavfiles/short1.wav"); /* WAVE�t�@�C������X�e���I�̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	pcm1_edge.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1_edge.bits = pcm0.bits; /* �ʎq�����x */
	pcm1_edge.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1_edge.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm1_edge.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	fe1 = 500.0 / (2 * pcm0.fs); /* �G�b�W���g��1 */
	fe2 = 2000.0 / (2 * pcm0.fs); /* �G�b�W���g��2 */
	delta = 900.0 / (2 * pcm0.fs); /* �J�ڑш敝 */

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

	b_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	b_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_adge = (double*)calloc(N, sizeof(double)); /* �������̊m�� */



	for (frame = 0; frame < number_of_frame; frame++)
	{
		offset = L * frame;

		/* x(n)��FFT */
		for (n = 0; n < N; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
		}
		for (n = 0; n < L/2; n++)
		{
			x_real[2 * n] = pcm0.sL[(offset/2) + n];
			x_real[(2 * n) + 1] = pcm0.sR[(offset / 2) + n];
		}
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

			y_adge[k] = 0;
		}

		IFFT(y_real, y_imag, N);

		/*�G�b�W��������*/
		//laplacian(y_real, y_adge, N);
		//edgestress(y_real, y_adge, N,enhancerate);


		/* �t�B���^�����O���ʂ̘A�� */
		for (n = 0; n < L * 2/2; n++)
		{
			if ((offset/2) + n < pcm1.length)
			{	
				pcm1.sL[(offset / 2) + n] += y_real[2 * n];
				pcm1.sR[(offset / 2) + n] += y_real[(2 * n) + 1];

				pcm1_edge.sL[(offset/2) + n] += y_adge[2 * n];			
				pcm1_edge.sR[(offset/2) + n] += y_adge[(2*n)+1];
			}
		}
	}

	/*�G�b�W��������*/
	laplacian(pcm1.sL, pcm1_edge.sL, pcm1_edge.length);
	laplacian(pcm1.sR, pcm1_edge.sR, pcm1_edge.length);
	//edgestress(y_real, y_adge, N,enhancerate);


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

	stereo_wave_write(&pcm1, "./wavfiles/out_short1.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */
	stereo_wave_write(&pcm1_edge, "./wavfiles/out_edge_short1.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */

	free(pcm0.sL); /* �������̉�� */
	free(pcm0.sR); /* �������̉�� */
	free(pcm1.sL); /* �������̉�� */
	free(pcm1.sR); /* �������̉�� */
	free(pcm1_edge.sL); /* �������̉�� */
	free(pcm1_edge.sR); /* �������̉�� */
	free(b); /* �������̉�� */
	free(w); /* �������̉�� */
	free(b_real); /* �������̉�� */
	free(b_imag); /* �������̉�� */
	free(x_real); /* �������̉�� */
	free(x_imag); /* �������̉�� */
	free(y_real); /* �������̉�� */
	free(y_imag); /* �������̉�� */
	free(y_adge); /* �������̉�� */

	//printf("counts =%d\n count = %d\n",counts,count);
	//printf("max = %lf\n ",max);
	printf("finish");

	return 0;
}
