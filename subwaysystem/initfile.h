#pragma once
typedef struct
{
	char* filename;
	char* directory;
	char* filetype;
} IN_FILE;

typedef struct
{
	char* orgfilename;
	char* directory;
	char* filetype;
	int fe1;
	int fe2;
	int delta;
	int done_bf_lap;
	int done_af_lap;
} OUT_FILE;

void init_infiles(IN_FILE *inputfile, char* directory, char* infilename, char* filetype) {
	//inputfile ‚ÌÝ’è
	inputfile->filename = infilename;
	inputfile->directory = directory;
	inputfile->filetype = filetype;
}

void init_outfiles(IN_FILE *inputfile, OUT_FILE *outputfiles, int fe1, int fe2, int delta) {
	//4‚Â‚Ì outputfile ‚ÌÝ’è
	for (int i = 0; i < 4; i++) {
		(outputfiles + i)->orgfilename = inputfile->filename;
		(outputfiles + i)->directory = inputfile->directory;
		(outputfiles + i)->filetype = inputfile->filetype;
		(outputfiles + i)->fe1 = fe1;
		(outputfiles + i)->fe2 = fe2;
		(outputfiles + i)->delta = delta;

		//{BP , BP.alap , blap.BP , blap.BP.alap}
		if (i < 2) {
			(outputfiles + i)->done_bf_lap = 0;
		}
		else {
			(outputfiles + i)->done_bf_lap = 1;
		}

		if (i % 2 == 0) {
			(outputfiles + i)->done_af_lap = 0;
		}
		else {
			(outputfiles + i)->done_af_lap = 1;
		}
	}
}
size_t getinfsize(IN_FILE *inputfile) {
	size_t infsize = strlen(inputfile->directory) + strlen(inputfile->filename) + strlen(inputfile->filetype) + 1;
	return infsize;
}

size_t getoutfsize(OUT_FILE *outputfile) {
	size_t outfsize = strlen(outputfile->directory) + 4 + strlen(outputfile->orgfilename) + outputfile->done_bf_lap * 6 + 18 + outputfile->done_af_lap * 6 + strlen(outputfile->filetype) + 1;
	return outfsize;
}

void getinfilename(IN_FILE *inputfile, char* filename) {

	strcpy(filename, inputfile->directory);
	strcat(filename, inputfile->filename);
	strcat(filename, inputfile->filetype);

}

void getoutfilename(OUT_FILE *outputfile, char* filename) {
	char frequency[19];
	sprintf(frequency, "_%d-%d(%d)", outputfile->fe1, outputfile->fe2, outputfile->delta);

	strcpy(filename, outputfile->directory);
	strcat(filename, "out_");
	strcat(filename, outputfile->orgfilename);
	strcat(filename, frequency);
	if (outputfile->done_bf_lap) {
		strcat(filename, "_bflap");
	}
	if (outputfile->done_af_lap) {
		strcat(filename, "_aflap");
	}
	strcat(filename, outputfile->filetype);
}