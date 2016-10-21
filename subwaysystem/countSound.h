#pragma once
#include "wave.h"

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