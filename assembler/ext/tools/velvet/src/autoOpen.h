#ifndef AUTOOPEN_H_
#define AUTOOPEN_H_

#include "globals.h"
#include "utility.h"

typedef struct {
	int pid;
	FILE* file;
	char const* decompressor;
	int first_char;
} AutoFile;

AutoFile* openFileAuto(char*filename);

void closeFileAuto(AutoFile* autoFile);

#endif
