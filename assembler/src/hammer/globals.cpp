//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "kmer_stat.hpp"
#include "globals.hpp"
#include "config_struct_hammer.hpp"

void Globals::writeBlob( const char * fname ) {
	std::ofstream ofs( fname );
	ofs << blob_max_size << "\n" << blob_size << "\n";
	ofs.write(blob, blob_size); ofs << "\n";
	ofs.write(blobquality, blob_size); ofs << "\n";
	ofs.close();
}


void Globals::readBlob( const char * fname ) {
	if (blob != NULL) delete [] blob;
	if (blobquality != NULL) delete [] blobquality;
	// if (blobhash != NULL) delete [] blobhash;

	FILE * f = fopen( fname, "r" );
	assert( fscanf(f, "%lu\n", &blob_max_size) != EOF );
	assert( fscanf(f, "%lu\n", &blob_size) != EOF );
	blob = new char[blob_max_size];
	blobquality = new char[blob_max_size];
	assert( fscanf(f, "%s\n", blob ) != EOF );
	assert( fscanf(f, "%s\n", blobquality ) != EOF );
	fclose(f);
}
