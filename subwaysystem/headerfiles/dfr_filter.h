#pragma once
//���ԗ̈�Ŏg�p��
void laplacian(double orgdatas[], double edgedatas[],int arraysize) {
	edgedatas[0] = 0;
	for (int i = 1; i < arraysize - 1; i++) {
		edgedatas[i] = orgdatas[i - 1] - 2 * orgdatas[i] + orgdatas[i + 1];
	}
	edgedatas[arraysize-1] = 0;
}

void normalize_array(double *orgarray,int datasize,double max) {
	int index;
	if (max != 0) {
		for (index = 0; index < datasize; index++) {
			orgarray[index] /= max;
		}
	}
	else {
		printf("error!!!!\nmax = 0\n");
	}
}


void edgestress(double *orgdatas, double *edgedatas, int arraysize,double rate) {
	for (int i = 0; i < arraysize; i++)
	{
		orgdatas[i] += rate*edgedatas[i];
	}
}

double getMax(double array[],int arraysize) {
	double max = 0,data = 0;
	for (int i = 0; i < arraysize ; i++) {
		data = fabs(array[i]);
		if (max < data) {
			max = data;
		//	printf("%lf\n", max);	//�f�o�b�O�R�[�h
		}
	}
	printf("array_max = %lf\n", max);	//�f�o�b�O�R�[�h
	return max;
}
double getMin(double *array, int arraysize) {
	double min = 1,data = 0;
	for (int i = 0; i < arraysize ; i++) {
		data = fabs(array[i]);
		if (min > data) {
			min = data;
			//printf("%lf\n", min);	//�f�o�b�O�R�[�h
		}
	}
	printf("array_min = %lf\n", min);	//�f�o�b�O�R�[�h
	return min;
}

double getbeforeMax(double *array, int arraysize) {
	double max = 0;
	double temp = 0;
	for (int i = 1; i < arraysize - 1; i++) {
		if (max < array[i]) {
			temp = max;
			max = array[i];
			printf("%lf\n", max);	//�f�o�b�O�R�[�h
		}
		else if (max < -array[i]) {
			temp = max;
			max = -array[i];
			printf("%-lf\n", max);	//�f�o�b�O�R�[�h
		}
	}
	printf("array_bfmax = %lf\n", temp);	//�f�o�b�O�R�[�h
	return temp;
}
