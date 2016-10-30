#pragma once
#include <stdio.h>

void make_datfile(char *filename, double data[], int datasize){

	FILE *fp;
	fp = fopen(filename, "w");

	for (int i = 0; i < datasize; i++) {
		fprintf(fp, "%d %lf\n",i, data[i]);
	}
	fprintf(fp,"\n\n");
	fclose(fp);

}

void postscript_datfile(char *filename, double data[],int datasize,int index){

	FILE *fp;
	fp = fopen(filename, "a");

	for (int i = 0; i < datasize; i++) {
		fprintf(fp, "%d %lf\n", index + i, data[i]);
	}	
	
	fclose(fp);	

}


//inputfile ‚Ìƒf[ƒ^”‚ðcompressionrate‚¾‚¯ˆ³k
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

