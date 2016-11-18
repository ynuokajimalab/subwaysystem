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

double getsec(int frame, int datanum, double fs);

int main(void)
{
	STEREO_PCM pcm0, pcm1;
	int m, n, J, k, N, L, offset, frame, number_of_frame, count, maxindex;
	double fe1, fe2, delta, *b, *b_real, *b_imag, *w, *x_real, *x_imag, *A, *y_real, *y_imag, max, average, sum, dispersion_ex, dispersion;

	stereo_wave_read(&pcm0, "short1.WAV"); /* WAVE�t�@�C������X�e���I�̉��f�[�^����͂��� */

	pcm1.fs = pcm0.fs; /* �W�{�����g�� */
	pcm1.bits = pcm0.bits; /* �ʎq�����x */
	pcm1.length = pcm0.length; /* ���f�[�^�̒��� */
	pcm1.sL = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */
	pcm1.sR = (double*)calloc(pcm1.length, sizeof(double)); /* �������̊m�� */

	fe1 = 1000.0 / (pcm0.fs * 2); /* �G�b�W���g�� */
	fe2 = 2000.0 / (pcm0.fs * 2); /* �G�b�W���g�� */
	delta = 800.0 / (pcm0.fs * 2); /* �J�ڑш敝 */

	J = (int)(3.1 / delta + 0.5) - 1; /* �x����̐� */
	if (J % 2 == 1)
	{
		J++; /* J+1����ɂȂ�悤�ɒ������� */
	}

	L = 256; /* �t���[���̒��� */
	N = 512; /* DFT�̃T�C�Y */

	number_of_frame = 2 * pcm0.length / L; /* �t���[���̐� */
	/*printf("�t���[�����F%d\n",number_of_frame);*/

	b = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */
	w = (double*)calloc((J + 1), sizeof(double)); /* �������̊m�� */

	Hanning_window(w, (J + 1)); /* �n�j���O�� */

	FIR_BPF(fe1, fe2, J, b, w); /* FIR�t�B���^�̐݌v */

	b_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	b_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	x_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	A = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_real = (double*)calloc(N, sizeof(double)); /* �������̊m�� */
	y_imag = (double*)calloc(N, sizeof(double)); /* �������̊m�� */

	count = 0;
	max = A[0];

	for (frame = 0; frame < number_of_frame; frame++)
	{
		sum = 0.0;
		dispersion_ex = 0.0;
		offset = L * frame;
		max = 0.0;
		maxindex = 0;

		/* x(n)��FFT */
		for (n = 0; n < N; n++)
		{
			x_real[n] = 0.0;
			x_imag[n] = 0.0;
		}
		for (n = 0; n < L / 2; n++)
		{
			x_real[2 * n] = pcm0.sL[(offset / 2) + n];
			x_real[2 * n + 1] = pcm0.sR[(offset / 2) + n];
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

		/* �U���X�y�N�g�� */
		for (k = 0; k < N; k++)
		{
			A[k] = sqrt(x_real[k] * x_real[k] + x_imag[k] * x_imag[k]);
		}

		/* �t�B���^�����O */
		for (k = 0; k < N; k++)
		{
			y_real[k] = x_real[k] * b_real[k] - x_imag[k] * b_imag[k];
			y_imag[k] = x_imag[k] * b_real[k] + x_real[k] * b_imag[k];
		}

		IFFT(y_real, y_imag, N);

		/* �t�B���^�����O���ʂ̘A�� */
		for (n = 0; n < L; n++)
		{
			if (offset / 2 + n < pcm1.length)
			{
				pcm1.sL[offset / 2 + n] += y_real[2 * n];
				pcm1.sR[offset / 2 + n] += y_real[2 * n + 1];
			}
		}

		for (k = 0; k < N; k++)
		{
			sum += A[k];
		}
		average = sum / N;
		for (k = 0; k < N; k++)
		{
			dispersion_ex += (average - A[k])*(average - A[k]);
		}
		dispersion = dispersion_ex / N;
		/*printf("dispersion=%lf\n", dispersion);
		printf("average=%lf\n",average);
*/
		for (k = 0; k < N; k++)
		{
			if (A[k] > 52) {
				count++;
				printf("frame=%d,A[%d]=%lf\n", frame, k, A[k]);
			}
		}

		/*for (k = 0; k < N; k++) {
			if (A[k] > max)
			{
				max = A[k];
				maxindex = k;
			}
		}
		if (A[maxindex] > 40)
		{
				count++;
				printf("frame = %d, time = %lf, max:A[%d]=%lf\n",frame,getsec(frame,L,pcm0.fs), maxindex, A[maxindex]);
		}*/

	}
	printf("�񐔁F%d\n", count);
	/*printf("fin_time = %lf\n",getsec(number_of_frame,L,pcm0.fs));*/

	stereo_wave_write(&pcm1, "try.WAV"); /* WAVE�t�@�C���ɃX�e���I�̉��f�[�^���o�͂��� */

	free(pcm0.sL); /* �������̉�� */
	free(pcm0.sR); /* �������̉�� */
	free(pcm1.sL); /* �������̉�� */
	free(pcm1.sR); /* �������̉�� */
	free(w); /* �������̉�� */
	free(x_real); /* �������̉�� */
	free(x_imag); /* �������̉�� */
	free(A); /* �������̉�� */
	free(y_real); /* �������̉�� */
	free(y_imag); /* �������̉�� */

	return 0;
}


double getsec(int frame, int datanum, double fs) {
	double sec;

	sec = ((double)(datanum*(frame + 1)) / (double)(4 * fs));

	return 0;
}