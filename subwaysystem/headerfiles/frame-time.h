#pragma once
int getframe(double second, int framedata, int SamplesPerSec) {
	int frame = ((int)((4 * SamplesPerSec)*second / framedata) - 1);
	if (frame < 0) {
		frame = 0;
	}
	return frame;
}

double getSecond(int offset, int SamplesPerSec) {
	return ((double)offset/(double)SamplesPerSec);
}

int getstereoindex(double time, int SamplesPerSec) {
	return((int)(time * SamplesPerSec*2));
}

void getpcm(STEREO_PCM *orgpcm,STEREO_PCM *newpcm, int index1, int index2) {

	int i;

	for (i = index1/2; i < index2/2; i++)
	{
		if (orgpcm->length < index2/2) {
			printf("指定された時刻でのデータはありません\n");
			return;
		}
		newpcm->sL[i - index1 /2] = orgpcm->sL[i];
		newpcm->sR[i - index1 /2] = orgpcm->sR[i];
	}


}