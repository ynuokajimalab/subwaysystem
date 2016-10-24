#pragma once
#include "wave.h"

int countsound(char *filename, double judge) {
	FILE *fp;
	int count, index;
	double data;

	count = 0;
	fp = fopen(filename, "r");

	while ((fscanf(fp, "%d %lf", &index, &data)) != EOF) {
		if (data >= judge || data <= -judge)
			count++;
	}

	fclose(fp);
	return count;
}

int judgeSounnd(double sound[] ,int N, double threshold) {
	int count = 0;
	for (int i = 0; i <N; i++)
	{
		if (sound[i] >= threshold || sound[i] <= -threshold){
			count = 1;
			break;
		}
	}
	return count;
}
