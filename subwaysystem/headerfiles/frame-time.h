#pragma once
int getframe(double second, int framedata, int SamplesPerSec) {
	int frame = ((int)((4 * SamplesPerSec)*second / framedata) - 1);
	if (frame < 0) {
		frame = 0;
	}
	return frame;
}

double getsecond(int frame, int framedata, int SamplesPerSec) {
	return ((double)(framedata *(frame + 1))/(double)(4 * SamplesPerSec));
}

int getstereoindex(double time, int SamplesPerSec) {
	return((int)(time * SamplesPerSec*2));
}

void getpcm(STEREO_PCM *orgpcm,STEREO_PCM *newpcm, double time1, double time2) {
	int index1, index2,i;
	index1 = getstereoindex(time1, orgpcm->fs);
	index2 = getstereoindex(time2, orgpcm->fs);

	for (i = index1/2; i < index2/2; i++)
	{
		newpcm->sL[i - index1 /2] = orgpcm->sL[i];
		newpcm->sR[i - index1 /2] = orgpcm->sR[i];
	}


}