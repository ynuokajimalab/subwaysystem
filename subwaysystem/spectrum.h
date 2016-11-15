#pragma once
//仮！！
void setnoise(double n_real[], double n_imag[], double x_real[], int size, int frame) {
	for (int k = 0; k < size; k++) {
		n_real[k] = (n_real[k] * frame + x_real[k]) / (frame + 1);
		n_imag[k] = 0.0;
	}
}

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


void subtruction(double org_amp[], double noize_amp[], int size) {
	for (int k = 0; k < size; k++) {
		org_amp[k] -= noize_amp[k];
		if (org_amp[k] < 0.0) {
			org_amp[k] = 0.0;
		}
	}
}

void do_subtruction(double s_real[], double s_imag[], double n_real[], double n_imag[], int size) {
	double *s_amp, *s_phas, *n_amp, *n_phas;
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
	free(s_amp);
	free(s_phas);
	free(n_amp);
	free(n_phas);
}