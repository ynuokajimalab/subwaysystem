#pragma warning(disable:4996)

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "./headerfiles/initfile.h"
#include "./headerfiles/wave.h"
#include "./headerfiles/window_function.h"
#include "./headerfiles/sinc.h"
#include "./headerfiles/fir_filter.h"
#include "./headerfiles/iir_filter.h"
#include "./headerfiles/fft.h"
#include "./headerfiles/countSound.h"
#include "./headerfiles/make_datfile.h"
#include "./headerfiles/dfr_filter.h"
#include "./headerfiles/spectrum.h"
#include "./headerfiles/frame-time.h"
#include "./headerfiles/statistics.h"

int main(void)
{
	STEREO_PCM org_pcm;
	char *infilename, *outfilename;
	FILE  *outfp;


	//入力ファイルのデータ
	char orgfile[] = "short2";
	char filetype[] = ".wav";
	char directory[] = "./wavfiles/";
	infilename = (char*)calloc(30, sizeof(char));
	outfilename = (char*)calloc(30, sizeof(char));

	
	strcpy(infilename, directory);
	strcat(infilename, orgfile);
	strcat(infilename,filetype);

	strcpy(outfilename, "");
	strcat(outfilename, orgfile);
	strcat(outfilename, ".txt");

	/* WAVEファイルからステレオの音データを入力する */
	stereo_wave_read(&org_pcm, infilename);


	outfp = fopen(outfilename, "w");

	for (int i = 0; i < org_pcm.length;i++) {
			fprintf(outfp, "%f\n", org_pcm.sL[i]);
			fprintf(outfp, "%f\n", org_pcm.sR[i]);
		}

	fclose(outfp);
	printf("ファイルを出力しました\n");
	printf("ファイル名：%s\n", outfilename);

	/* メモリの解放 */
	free(infilename);
	free(outfilename);

	printf("finish\n");

	return 0;
}
