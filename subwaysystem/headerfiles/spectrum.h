#pragma once
//âºÅIÅI

void getampphase(double real[], double image[], double amplitude[], double phase[], int size) {
	for (int k = 0; k < size; k++) {
		amplitude[k] = sqrt(real[k] * real[k] + image[k] * image[k]);
		if (image[k] != 0.0 && real[k] != 0.0) {
			phase[k] = atan2(image[k], real[k]);
		}
	}
}

void getrealimage(double real[], double image[], double amplitude[], double phase[], int size) {
	for (int k = 0; k < size; k++) {
		real[k] = amplitude[k] * cos(phase[k]);
		image[k] = amplitude[k] * sin(phase[k]);
	}
}

void setnoise(double n_real[], double n_imag[],double n_amp[], double n_phas[], double w[],STEREO_PCM pcm0, int FFTsize,int framelength, int noise_frame1, int noise_frame2) {
	int k,n,offset;
	int frame;

	for (frame = noise_frame1; frame < noise_frame2;frame++) {
		offset = (framelength / 2) * frame;
		for (n = 0; n < FFTsize / 2; n++)
		{
			n_real[2 * n] += pcm0.sL[(offset / 2) + n];
			n_real[(2 * n) + 1] += pcm0.sR[(offset / 2) + n];
		}
	}

	for (k = 0; k < FFTsize; k++) {
		n_real[k] = (n_real[k] * w[k])/(noise_frame2-noise_frame1+1);
		//n_real[k] = (n_real[k] * w[k]) / sqrt(noise_frame2 - noise_frame1 + 1);

		n_imag[k] = 0.0;
	}
	FFT(n_real, n_imag, FFTsize);
	getampphase(n_real,n_imag,n_amp,n_phas,FFTsize);
	IFFT(n_real, n_imag, FFTsize);

	//for (k = 0; k < FFTsize; k++) {
	//	printf("n_amp[%d] = %lf\n",k, n_amp[k]);
	//}
}

void subtruction(double org_amp[], double noize_amp[], int size) {
	for (int k = 0; k < size; k++) {
		//printf("%d::%lf - %lf = ",k,org_amp[k],noize_amp[k]);
		//org_amp[k] -= noize_amp[k];
		org_amp[k] -= 3.0;//test cord
		if (org_amp[k] < 0.0) {
			org_amp[k] = 0.0;
		}
		//printf("%lf\n",org_amp[k]);
	}
}

//éûä‘óÃàÊÇ≈ì¸óÕ
void do_subtruction(double s_real[], double s_imag[], double n_real[], double n_imag[],double w[], int size) {
	double *s_amp, *s_phas, *n_amp, *n_phas;

	for (int i = 0; i < size;i++) {
		s_real[i] *= w[i];
	}

	FFT(n_real, n_imag, size);
	FFT(s_real, s_imag, size);

	s_amp = (double*)calloc(size, sizeof(double));
	s_phas = (double*)calloc(size, sizeof(double));
	n_amp = (double*)calloc(size, sizeof(double));
	n_phas = (double*)calloc(size, sizeof(double));
	for (int k = 0; k < size; k++) {
		s_amp[k] = 0.0;
		s_phas[k] = 0.0;
		n_amp[k] = 0.0;
		n_phas[k] = 0.0;
	}
	getampphase(s_real, s_imag, s_amp, s_phas, size);
	getampphase(n_real, n_imag, n_amp, n_phas, size);
	subtruction(s_amp, n_amp, size);
	getrealimage(s_real, s_imag, s_amp, s_phas, size);

	IFFT(n_real, n_imag, size);
	IFFT(s_real, s_imag, size);

	free(s_amp);
	free(s_phas);
	free(n_amp);
	free(n_phas);
}
