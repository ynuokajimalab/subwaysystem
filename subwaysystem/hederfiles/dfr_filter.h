#pragma once
//時間領域で使用可
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

double getMax(double *array,int arraysize) {
	double max = 0;
	for (int i = 1; i < arraysize - 1; i++) {
		if (max < array[i]) {
			max = array[i];
			printf("%lf\n", max);	//デバッグコード
		}else if (max < -array[i]) {
			max = -array[i];
			printf("%lf\n",max);	//デバッグコード
		}
	}
	printf("array_max = %lf\n", max);	//デバッグコード
	return max;
}

