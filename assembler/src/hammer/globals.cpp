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
	// if (blobhash != NULL) delete [] blobhash;

	FILE * f = fopen( fname, "r" );
	assert( fscanf(f, "%lu\n", &blob_max_size) != EOF );
	assert( fscanf(f, "%lu\n", &blob_size) != EOF );
	blob = new char[blob_max_size];
	blobquality = new char[blob_max_size];
	assert( fscanf(f, "%s\n", blob ) != EOF );
	assert( fscanf(f, "%s\n", blobquality ) != EOF );
	fclose(f);

	// precompute hashes
	//Globals::blobhash = new uint64_t[ Globals::blob_max_size ];
	//KMerNo::precomputeHashes();
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

void Globals::readKMerFiles( const char * f_num, const char * f_solid, const char * f_bad, std::vector<KMerCount*> * kmers ) {
	// get number of kmers
	ifstream ifs(f_num);
	string num; getline(ifs, num); hint_t totalKMers;
	assert(sscanf(num.c_str(), "%lu", &totalKMers) != EOF);
	ifs.close();

	for ( size_t i=0; i < kmers->size(); ++i ) delete (*kmers)[i]; kmers->clear();
	kmers->resize(totalKMers, NULL);
	ifs.open(f_solid);
	unsigned long int start; unsigned int count; unsigned long int index; unsigned long int changeto; double qual; double clust; string kmer; string line;
	while (ifs.good()) {
		getline(ifs, kmer);
		if (kmer.size() != K) break; // we're at the end of file
		getline(ifs, line);
		if (sscanf(line.c_str(), ">%lu\tgood singleton\tind=%lu\tcnt=%u\ttql=%lf", &start, &index, &count, &qual) != EOF ) {
			cout << "good! " << line << endl;
			cout << start << " " << index << " " << count << " " << qual << " " << (qual > Globals::iterative_reconstruction_threshold) << endl;
			(*kmers)[index] = new KMerCount( PositionKMer(start),
					  KMerStat(count, (qual > Globals::iterative_reconstruction_threshold ? KMERSTAT_GOODITER : KMERSTAT_GOOD), qual) );
		} else if (sscanf(line.c_str(), ">%lu\tcenter clust=%lf\tind=%lu\tcnt=%u\ttql=%lf", &start, &clust, &index, &count, &qual) != EOF) {
			cout << "center! " << line << endl;
			cout << start << " " << index << " " << clust << " " << count << " " << qual << " " << (clust > Globals::iterative_reconstruction_threshold) << endl;
			(*kmers)[index] = new KMerCount( PositionKMer(start),
					  KMerStat(count, (clust > Globals::iterative_reconstruction_threshold ? KMERSTAT_GOODITER : KMERSTAT_GOOD), qual) );
		} else {
			assert(false);
		}
	}
	ifs.close();
	ifs.open(f_bad);
	while (ifs.good()) {
		getline(ifs, kmer);
		if (kmer.size() != K) break; // we're at the end of file
		getline(ifs, line);
		if (sscanf(line.c_str(), ">%lu\tbad singleton\tind=%lu\tcnt=%u\ttql=%lf", &start, &index, &count, &qual) != EOF ) {
			(*kmers)[index] = new KMerCount( PositionKMer(start),
					  KMerStat(count, KMERSTAT_BAD, qual) );
		} else if (sscanf(line.c_str(), ">%lu\tcenter of bad clust=%lf\tind=%lu\tcnt=%u\ttql=%lf", &start, &clust, &index, &count, &qual) != EOF) {
			cout << "bad center! " << line << endl;
			cout << start << " " << index << " " << clust << " " << count << " " << qual << " " << (clust > Globals::iterative_reconstruction_threshold) << endl;
			(*kmers)[index] = new KMerCount( PositionKMer(start),
					  KMerStat(count, KMERSTAT_BAD, qual) );
		} else if (sscanf(line.c_str(), ">%lu\tpart of cluster %lu\tclust=%lf\tind=%lu\tcnt=%u\ttql=%lf", &start, &changeto, &clust, &index, &count, &qual) != EOF) {
			(*kmers)[index] = new KMerCount( PositionKMer(start),
					  KMerStat(count, changeto, qual) );
		} else {
			assert(false);
		}
	}
	ifs.close();
}


void Globals::writeKMerHashMap( const char * fname, const KMerNoHashMap & hm ) {
	std::ofstream ofs( fname );
	for( KMerNoHashMap::const_iterator it = hm.begin(); it != hm.end(); ++it ) {
		ofs << it->first.index << "\t" << it->second->first.start() << "\t" << it->second->second.count << "\t" << it->second->second.totalQual << "\t" << it->second->second.changeto << "\n";
	}
	ofs.close();
}

void Globals::readKMerHashMap( const char * fname, KMerNoHashMap * hm, std::vector<KMerCount*> * kmers ) {
	for ( size_t i=0; i < kmers->size(); ++i ) delete (*kmers)[i]; kmers->clear(); hm->clear();
	ifstream ifs(fname);
	hint_t index; uint32_t count; double totalQual; hint_t changeto; hint_t start; std::string line;
	while (ifs.good()) {
		getline(ifs, line);
		if (sscanf(line.c_str(), "%lu\t%lu\t%u\t%lf\t%lu", &index, &start, &count, &totalQual, &changeto) != EOF ) {
			KMerCount * kmc = new KMerCount( PositionKMer(start), KMerStat(count, changeto, totalQual) );
			hm->insert(make_pair( KMerNo(index, 1), kmc) );
			kmers->push_back(kmc);
		}
	}
	ifs.close();
}

