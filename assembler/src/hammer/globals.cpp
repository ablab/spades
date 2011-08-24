#include "kmer_stat.hpp"
#include "position_kmer.hpp"
#include "globals.hpp"

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
	if (blobhash != NULL) delete [] blobhash;

	FILE * f = fopen( fname, "r" );
	assert( fscanf(f, "%lu\n", &blob_max_size) != EOF );
	assert( fscanf(f, "%lu\n", &blob_size) != EOF );
	blob = new char[blob_max_size];
	blobquality = new char[blob_max_size];
	assert( fscanf(f, "%s\n", blob ) != EOF );
	assert( fscanf(f, "%s\n", blobquality ) != EOF );
	fclose(f);

	// precompute hashes
	Globals::blobhash = new uint64_t[ Globals::blob_max_size ];
	KMerNo::precomputeHashes();
}

void Globals::writeKMerCounts( const char * fname, const std::vector<KMerCount*> & kmers ) {
	std::ofstream ofs( fname );
	for ( size_t i=0; i < kmers.size(); ++i ) {
		ofs << kmers[i]->first.str() << "\t" << kmers[i]->first.start() << "\t" << kmers[i]->second.count << "\t" << kmers[i]->second.totalQual << "\n";
	}
	ofs.close();
}

void Globals::readKMerCounts( const char * fname, std::vector<KMerCount*> * kmers ) {
	kmers->clear();
	FILE * f = fopen( fname, "r" );
	unsigned long int start; unsigned int count; unsigned long int startlast = -1; double qual; char tmp[K+10];
	while (!feof(f)) {
		assert( fscanf(f, "%s\t%lu\t%u\t%lf", tmp, &start, &count, &qual) != EOF );
		if (start != startlast) {
			kmers->push_back( new KMerCount( PositionKMer(start), KMerStat(count, KMERSTAT_GOOD, qual) ) );
			startlast = start;
		}
	}
	fclose(f);
}

