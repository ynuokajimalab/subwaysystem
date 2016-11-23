#pragma once
typedef struct
{
	char* filename;
	char* directory;
	char* filetype;
} FILENAME;

void init_file(FILENAME *file, char* directory, char* infilename, char* filetype) {
	//file ‚ÌÝ’è
	strcpy(file->filename,infilename);
	file->directory = directory;
	file->filetype = filetype;
}

size_t getfnsize(FILENAME file) {
	size_t infsize = strlen(file.directory) + strlen(file.filename) + strlen(file.filetype) + 1;
	return infsize;
}


void getfilename(FILENAME file, char* filename) {

	int i = 0;
	strcpy(filename, file.directory);
	strcat(filename, file.filename);
	strcat(filename, file.filetype);
}