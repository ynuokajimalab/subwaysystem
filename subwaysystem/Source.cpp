#pragma warning(disable:4996)

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "wave.h"
#include "window_function.h"
#include "sinc.h"
#include "fir_filter.h"
#include "fft.h"
#include "countSound.h"
#include "make_datfile.h"


int main(void)
{
	STEREO_PCM pcm0, pcm1;
	int n, m, k, J, L, N, offset, frame, number_of_frame,countsound,compressionrate,filehasdata;
	double fe1,fe2, delta, threshold, *b, *w, *b_real, *b_imag, *x_real, *x_imag, *y_real, *y_imag;


	stereo_wave_read(&pcm0, "short1.wav"); /* WAVE�t�@�C������X�e���I�̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	fe1 = 150.0 / pcm0.fs; /* �G�b�W���g��1 */
	fe2 = 450.0 / pcm0.fs; /* �G�b�W���g��2 */
	delta = 450.0 / pcm0.fs; /* �J�ڑш敝 */

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

	countsound = 0; /* ����臒l�𒴂�����*/
	threshold = 10.0;/* 臒l*/

	filehasdata = 0; /*�t�@�C����V������邩���肷��t���O*/
	compressionrate = 500; /*�f�[�^�̈��k��(-�{)*/


	number_of_frame = pcm0.length*2 / L; /* �t���[���̐� */


	b_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	b_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */

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
		}

		countsound += judgeSounnd(y_real,N,threshold);
		IFFT(y_real, y_imag, N);

		if (filehasdata != 0) {
			postscript_datfile("allsound.dat", y_real, N,offset);
			postscript_datfile("orgsound.dat", x_real, N, offset);

		}
		else {
			make_datfile("allsound.dat",y_real,N);
			make_datfile("orgsound.dat", x_real, N);
			filehasdata = 1;
		}

		if (frame == number_of_frame/2)
			make_datfile("onepointsound.dat",y_real,N);


		/* �t�B���^�����O���ʂ̘A�� */
		for (n = 0; n < L * 2/2; n++)
		{
			if ((offset/2) + n < pcm1.length)
			{	
				pcm1.sL[(offset/2) + n] += y_real[2 * n];			
				pcm1.sR[(offset/2) + n] += y_real[(2*n)+1];
			}
		}
	}
	double max1 = sampling_data("allsound.dat", "samplingsound.dat", compressionrate);
	normalize_data("samplingsound.dat","nomalizesound.dat",max1);

	stereo_wave_write(&pcm1, "out_short1_BPFver.wav"); /* WAVE�t�@�C���Ƀ��m�����̉��f�[�^���o�͂��� */

	free(pcm0.sL); /* �������̉�� */
	free(pcm0.sR); /* �������̉�� */
	free(pcm1.sL); /* �������̉�� */
	free(pcm1.sR); /* �������̉�� */
	free(b); /* �������̉�� */
	free(w); /* �������̉�� */
	free(b_real); /* �������̉�� */
	free(b_imag); /* �������̉�� */
	free(x_real); /* �������̉�� */
	free(x_imag); /* �������̉�� */
	free(y_real); /* �������̉�� */
	free(y_imag); /* �������̉�� */

	printf("countsound =%d\n",countsound);
	printf("finish");

	return 0;
}
