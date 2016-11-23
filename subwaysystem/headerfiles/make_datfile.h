#pragma once
#include <stdio.h>

typedef struct
{
	double time;
	int freequency;
	double power;
} SPECTROGRAM;


void make_datfile(char *filename){

	FILE *fp;
	fp = fopen(filename, "w");
	printf("%s:データファイル作成\n",filename);
	fclose(fp);

}

void postscript_datfile(char *filename,SPECTROGRAM spectorogram[], int samplenumber){

	FILE *fp;
	fp = fopen(filename, "a");

	for (int i = 0; i < samplenumber; i++) {
		fprintf(fp, "%lf %d %lf\n",spectorogram[i].time,spectorogram[i].freequency,spectorogram[i].power);
	}	
	
	fclose(fp);	

}

void makehistgram(char *filename,double histgram[],double dn, int histgramsize) {

	FILE *fp;
	fp = fopen(filename, "w");

	for (int i = 0; i < histgramsize; i++) {
		fprintf(fp, "%lf %lf\n", i*dn, histgram[i]);
	}

	fclose(fp);

}



//inputfile のデータ数をcompressionrateだけ圧縮
double sampling_data(char *inputfilename,char *outputfilename,int compressionrate){
	int index;
	double data;
	FILE *infp,*outfp;
	double max;

	infp = fopen(inputfilename, "r");
	outfp = fopen(outputfilename, "w");
	index = 0;
	max = 0;

		
	while((fscanf(infp, "%d %lf", &index, &data)) != EOF){

		if (index%compressionrate == 0) {
			if (max < data) {
				max = data;
			}
			fprintf(outfp,"%d %f\n",index,data);
		}
	}

	fclose(infp);
	fclose(outfp);

	return max;

}

void normalize_data(char *inputfilename, char *outputfilename, double max) {
	double data;
	int index;
	FILE *infp, *outfp;

	infp = fopen(inputfilename, "r");
	outfp = fopen(outputfilename, "w");


	while ((fscanf(infp, "%d %lf", &index, &data)) != EOF) {
		data /= max;
		fprintf(outfp, "%d %lf\n", index, data);
	}

	fclose(infp);
	fclose(outfp);

}

