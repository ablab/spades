#include "position_kmer.hpp"

void PositionKMer::writeBlob( const char * fname ) {
	ofstream ofs( fname );
	ofs << blob_max_size << "\n" << blob_size << "\n";
	ofs.write(blob, blob_size); ofs << "\n";
	ofs.write(blobquality, blob_size); ofs << "\n";
	for (hint_t i=0; i < blob_size; ++i) {
		ofs << blobkmers[i] << "\n";
	}
	ofs.close();
}


void PositionKMer::readBlob( const char * fname ) {
	if (blob != NULL) delete [] blob;
	if (blobquality != NULL) delete [] blobquality;
	if (blobkmers != NULL) delete [] blobkmers;
	if (blobprob != NULL) delete [] blobprob;

	FILE * f = fopen( fname, "r" );
	assert( fscanf(f, "%lu\n", &blob_max_size) != EOF );
	assert( fscanf(f, "%lu\n", &blob_size) != EOF );
	blob = new char[blob_max_size];
	blobquality = new char[blob_max_size];
	blobkmers = new hint_t[blob_max_size];
	assert( fscanf(f, "%s\n", blob ) != EOF );
	assert( fscanf(f, "%s\n", blobquality ) != EOF );
	hint_t tmp;
	for (hint_t i=0; i < blob_size; ++i) {
		if (feof(f)) { cout << "Not enough blobkmers!" << endl; break; }
		assert( fscanf(f, "%lu\n", &tmp) != EOF );
		blobkmers[i] = tmp;
	}
	fclose(f);
}

void PositionKMer::writeKMerCounts( const char * fname, const vector<KMerCount> & kmers ) {
	ofstream ofs( fname );
	for ( size_t i=0; i < kmers.size(); ++i ) {
		ofs << kmers[i].first.str() << "\t" << kmers[i].first.start() << "\t" << kmers[i].second.count << "\t" << kmers[i].second.totalQual << "\n";
	}
	ofs.close();
}

void PositionKMer::readKMerCounts( const char * fname, vector<KMerCount> * kmers ) {
	kmers->clear();
	FILE * f = fopen( fname, "r" );
	unsigned long int start; unsigned int count; unsigned long int startlast = -1; double qual; char tmp[K+10];
	while (!feof(f)) {
		assert( fscanf(f, "%s\t%lu\t%u\t%lf", tmp, &start, &count, &qual) != EOF );
		if (start != startlast) {
			kmers->push_back( make_pair( PositionKMer(start), KMerStat(count, KMERSTAT_GOOD, qual) ) );
			startlast = start;
		}
	}
	fclose(f);
}

