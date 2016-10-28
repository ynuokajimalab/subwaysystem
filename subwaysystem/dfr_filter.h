#pragma once
//ŽžŠÔ—Ìˆæ‚ÅŽg—p‰Â
void laplacian(double orgdatas[], double edgedatas[],int arraysize) {
	edgedatas[0] = 0;
	for (int i = 1; i < arraysize - 1; i++) {
		edgedatas[i] = orgdatas[i - 1] - 2 * orgdatas[i] + orgdatas[i + 1];
	}
	edgedatas[arraysize-1] = 0;
}

void edgestress(double *orgdatas, double *edgedatas, int arraysize,double rate) {
	for (int i = 0; i < arraysize; i++)
	{
		orgdatas[i] += rate*edgedatas[i];
	}
}