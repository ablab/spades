/*
 * hammer_tools.cpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#include "standard.hpp"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/format.hpp>
#include <time.h>
#include <sys/resource.h>
#include <iomanip>

#include "read/ireadstream.hpp"
#include "defs.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
#include "position_kmer.hpp"
#include "globals.hpp"
#include "config_struct_hammer.hpp"
#include "subkmers.hpp"
#include "hammer_tools.hpp"

// forking
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

string encode3toabyte (const string & s)  {
	string retval;
	char c = 48;
	int weight = 16;
	size_t i;
	for (i = 0; i < s.length(); i += 1) {
		if (i % 3 == 0) {
			c= 48;
			weight = 16;
		}
		c += weight * nt2num(s[i]);
		weight /= 4;
		if (i % 3 == 2) retval += c;
	}
	if (i % 3 != 0) retval += c;
	return retval;
}

void join_maps(KMerStatMap & v1, const KMerStatMap & v2) {
	KMerStatMap::iterator itf;
	for (KMerStatMap::const_iterator it = v2.begin(); it != v2.end(); ++it) {
		itf = v1.find(it->first);
		if (itf != v1.end()) {
			itf->second.count = itf->second.count + it->second.count;
			itf->second.totalQual = itf->second.totalQual * it->second.totalQual;
		} else {
			PositionKMer pkm = it->first;
			KMerStat kms( it->second.count, KMERSTAT_GOODITER, it->second.totalQual );
			v1.insert( make_pair( pkm, kms ) );
		}
	}
}

void print_time() {
	time_t rawtime;
	tm * ptm;
	time ( &rawtime );
	ptm = gmtime( &rawtime );
	cout << setfill('0') << "[ " << setw(2) << ptm->tm_hour << ":" << setw(2) << ptm->tm_min << ":" << setw(2) << ptm->tm_sec << " ] ";
}

void print_mem_usage() {
	static size_t pid = getpid();
	std::string str = (boost::format("pmap -d %d | grep writeable/private") % pid).str();
	std::cout << "==== MEM USAGE ==== " << std::endl;
	if (system(str.c_str())) { std::cout << "  System error!" << endl; }
}

void print_stats() {
	std::cout << "[";
	time_t rawtime;
	tm * ptm2;
	time ( &rawtime );
	ptm2 = gmtime( &rawtime );
	cout << " " << setfill('0') << setw(2) << ptm2->tm_hour << ":" << setw(2) << ptm2->tm_min << ":" << setw(2) << ptm2->tm_sec;

	rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	tm * ptm = gmtime( &ru.ru_utime.tv_sec );
	std::cout << " " << setw(2) << setfill('0')  << ptm->tm_hour << ":" << setw(2) << setfill('0') << ptm->tm_min << ":" << setw(2) << setfill('0') << ptm->tm_sec;
	if ( ru.ru_maxrss < 1024 * 1024 ) {
		std::cout << setw(4) << setfill(' ') << (ru.ru_maxrss / 1024) << "M ";
	} else {
		std::cout << setw(4) << setprecision(1) << fixed << setfill(' ') << (ru.ru_maxrss / (1024.0 * 1024.0) ) << "G ";
	}
	std::cout << "] ";
}


void HammerTools::ChangeNtoAinReadFiles() {
	std::vector<pid_t> pids;
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		string cur_filename = getFilename(cfg::get().input_working_dir, "reads.input", iFile);
		if (cfg::get().count_do) {
			pid_t cur_pid = vfork();
			if (cur_pid == 0) {
				TIMEDLN("  [" << getpid() << "] Child process for substituting Ns in " << Globals::input_filenames[iFile] << " starting.");
				string cmd = string("sed \'n;s/\\([ACGT]\\)N\\([ACGT]\\)/\\1A\\2/g;n;n\' ") + Globals::input_filenames[iFile].c_str() + " > " + cur_filename.c_str();
				int exitcode = system( cmd.c_str() );
				if (exitcode != 0) {
					TIMEDLN("  [" << getpid() << "] ERROR: finished with non-zero exit code " << exitcode);
				}
				_exit(0);
			}
			pids.push_back(cur_pid);
		}
		Globals::input_filenames[iFile] = cur_filename;
	}
	if (cfg::get().count_do) {
		int childExitStatus;
		for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
			waitpid(pids[iFile], &childExitStatus, 0);
		}
	}
}

void HammerTools::DecompressIfNeeded() {
	struct stat st;
	char f_cont[2];
	vector<pid_t> pIDs(Globals::input_filenames.size(), 0);
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		stat(Globals::input_filenames[iFile].c_str(), &st);
		std::ifstream ifs(Globals::input_filenames[iFile], std::ios::in | std::ios::binary );
		ifs.read(f_cont, 2);
		ifs.close();
		if ( (f_cont[0] == (char)0x1f) && (f_cont[1] == (char)0x8b) ) {
			string newFilename = HammerTools::getFilename(cfg::get().input_working_dir, "reads.input", iFile);
			string oldFilename = Globals::input_filenames[iFile];
			Globals::input_filenames[iFile] = newFilename;
			pIDs[iFile] = vfork();
			if (pIDs[iFile] == 0) {
				string systemcall = string("gunzip -c ") + oldFilename + string(" > ") + newFilename;
				TIMEDLN("  [" << getpid() << "] " << systemcall);
				if (system(systemcall.c_str())) {
					TIMEDLN("  [" << getpid() << "] System error with unzipping input files!");
				}
			_exit(0);
			}
		}
	}
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		if (pIDs[iFile] != 0) {
			int childExitStatus;
			waitpid(pIDs[iFile], &childExitStatus, 0);
		}
	}
}

hint_t HammerTools::EstimateTotalReadSize() {
	struct stat st;
	hint_t totalReadSize = 0;
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		stat(Globals::input_filenames[iFile].c_str(), &st);
		totalReadSize += st.st_size;
	}
	totalReadSize = totalReadSize / (2.5);
	return totalReadSize;
}


void HammerTools::InitializeSubKMerPositions() {
	ostringstream log_sstream;
	log_sstream.str("");
	Globals::subKMerPositions = new std::vector<uint32_t>(cfg::get().general_tau + 2);
	for (uint32_t i=0; i < (uint32_t)(cfg::get().general_tau + 1); ++i) {
		Globals::subKMerPositions->at(i) = (i * K / (cfg::get().general_tau + 1) );
		log_sstream << Globals::subKMerPositions->at(i) << " ";
	}
	Globals::subKMerPositions->at(cfg::get().general_tau + 1) = K;
	TIMEDLN("Hamming graph threshold tau=" << cfg::get().general_tau << ", k=" << K << ", subkmer positions = [ " << log_sstream.str() << "]" );
}

void HammerTools::ReadFileIntoBlob(const string & readsFilename, hint_t & curpos, hint_t & cur_read, bool reverse_complement) {
	TIMEDLN("Reading input file " << readsFilename);
	ireadstream irs(readsFilename, cfg::get().input_qvoffset);
	Read r;
	while (irs.is_open() && !irs.eof()) {
		irs >> r;
		size_t read_size = r.trimNsAndBadQuality(cfg::get().input_trim_quality);
		//cout << r.getSequenceString().c_str() << endl;
		if (read_size < K) continue;
		if ( reverse_complement ) r = !r;
		PositionRead pread(curpos, read_size, cur_read, false);
		Globals::pr->push_back(pread);
		for (uint32_t j = 0; j < read_size; ++j) {
			Globals::blob[curpos + j] = r.getSequenceString()[j];
			Globals::blobquality[curpos + j] = (char) (cfg::get().input_qvoffset + r.getQualityString()[j]);
		}
		//cout << "  read " << cur_read << "." << reverse_complement << ": pos=" << curpos << " size=" << read_size << endl;
		curpos += read_size;
		++cur_read;
	}
	irs.close();
}

void HammerTools::ReadAllFilesIntoBlob() {
	if (Globals::pr) Globals::pr->clear(); else Globals::pr = new vector<PositionRead>();
	hint_t curpos = 0;
	hint_t cur_read = 0;
	Globals::input_file_blob_positions.clear();
	Globals::input_file_blob_positions.push_back(0);
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		ReadFileIntoBlob(Globals::input_filenames[iFile], curpos, cur_read, false);
		Globals::input_file_blob_positions.push_back(cur_read);
	}
	Globals::revNo = cur_read;
	Globals::revPos = curpos;
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		ReadFileIntoBlob(Globals::input_filenames[iFile], curpos, cur_read, true);
	}
	Globals::blob_size = curpos;
	//std::cout << "  real blob size = " << curpos << endl;
	//std::cout << string(Globals::blob, curpos) << endl;
	//std::cout << "blob with revcomp:\n" << string(Globals::blob + Globals::revPos, curpos - Globals::revPos) << endl;
}

void HammerTools::CountAndSplitKMers(bool writeFiles) {
	TIMEDLN("Splitting kmer instances into files.");
	vector<boost::shared_ptr<FOStream> > ostreams(cfg::get().count_numfiles);
	if (writeFiles) {
		for (int i = 0; i < cfg::get().count_numfiles; ++i) {
			ostreams[i] = FOStream::init(getFilename( cfg::get().input_working_dir, Globals::iteration_no, "tmp.kmers", i ));
		}
	}
	Seq<K>::hash hash_function;
	//Seq<K>::less2 cmp_kmers;
	for(size_t i=0; i < Globals::revNo; ++i) {
		string s(Globals::blob        + Globals::pr->at(i).start(), Globals::pr->at(i).size());
		string q(Globals::blobquality + Globals::pr->at(i).start(), Globals::pr->at(i).size());
		for ( size_t j=0; j < Globals::pr->at(i).size(); ++j) q[j] = (char)(q[j] - cfg::get().input_qvoffset);
		ValidKMerGenerator<K> gen(s, q);
		while (gen.HasMore()) {
			//Seq<K> cur_kmer = gen.kmer();
			//Seq<K> cur_rckmer = !cur_kmer;
			//hint_t cur_pos = Globals::pr->at(i).start() + gen.pos() - 1;
			/*if (cmp_kmers(cur_rckmer, cur_kmer)) {
				cur_kmer = cur_rckmer;
				cur_pos = Globals::pr->at(i).getRCPosition(gen.pos());
				Globals::pr->at(i).setRCBit(gen.pos() - 1);
			}*/
			if (writeFiles) {
				boost::iostreams::filtering_ostream & cur_gz = ostreams[hash_function(cur_kmer) % cfg::get().count_numfiles]->fs;
				double correct_probability = 1 - gen.correct_probability();
				cur_gz << cur_pos << "\t" << correct_probability << "\t" << "\n";
			}
			gen.Next();
		}
	}
	ostreams.clear();
}


void HammerTools::CountKMersBySplitAndMerge() {
	CountAndSplitKMers(true);
	int count_num_threads = min( cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads );
	TIMEDLN("Kmer instances split. Starting merge in " << count_num_threads << " threads.");
	boost::shared_ptr<FOStream> kmerno_gz = FOStream::init(getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total"));
	hint_t kmer_num = 0;

	int merge_nthreads = min( cfg::get().general_max_nthreads, cfg::get().count_merge_nthreads);

	#pragma omp parallel for shared(kmerno_gz, kmer_num) num_threads(merge_nthreads)
	for ( int iFile=0; iFile < cfg::get().count_numfiles; ++iFile ) {
		KMerNoHashMap khashmap;
		boost::shared_ptr<FIStream> inStream = FIStream::init(getFilename( cfg::get().input_working_dir, Globals::iteration_no, "tmp.kmers", iFile ));
		ProcessKmerHashFile( inStream->fs, khashmap );
		#pragma omp critical
		{
			PrintProcessedKmerHashFile( kmerno_gz->fs, kmer_num, khashmap );
			TIMEDLN("Processed file " << iFile << " with thread " << omp_get_thread_num());
		}
	}

	TIMEDLN("Merge done. There are " << kmer_num << " kmers in total.");
}

hint_t HammerTools::IterativeExpansionStep(int expand_iter_no, int nthreads, const vector<KMerCount*> & kmers) {
	hint_t res = 0;

	// cycle over the reads, looking for reads completely covered by solid k-mers
	// and adding new solid k-mers on the fly
	#pragma omp parallel for shared(res) num_threads(nthreads)
	for (hint_t readno = 0; readno < Globals::revNo; ++readno) {
		// maybe this read has already been covered by solid k-mers
		if (Globals::pr->at(readno).isDone()) continue;

		const PositionRead & pr = Globals::pr->at(readno);
		const uint32_t read_size = pr.size();
		string seq = string(Globals::blob + pr.start(), read_size);
		vector< bool > covered_by_solid( read_size, false );
		vector< hint_t > kmer_indices( read_size, -1 );

		pair<int, hint_t > it = make_pair( -1, BLOBKMER_UNDEFINED );

		while ( (it = pr.nextKMerNo(it.first)).first > -1 ) {
			kmer_indices[it.first] = it.second;
			if ( kmers[it.second]->second.isGoodForIterative() ) {
				for ( size_t j = it.first; j < it.first + K; ++j )
					covered_by_solid[j] = true;
			}
		}

		bool isGood = true;
		for ( size_t j = 0; j < read_size; ++j ) {
			if ( !covered_by_solid[j] ) { isGood = false; break; }
		}
		if ( !isGood ) continue;

		// ok, now we're sure that everything is covered
		// first, set this read as already done
		Globals::pr->at(readno).done();

		// second, mark all k-mers as solid
		for (size_t j = 0; j < read_size; ++j) {
			if ( kmer_indices[j] == (hint_t)-1 ) continue;
			if ( !kmers[kmer_indices[j]]->second.isGoodForIterative() &&
				 !kmers[kmer_indices[j]]->second.isMarkedGoodForIterative() ) {
				#pragma omp critical
				{
				++res;
				kmers[kmer_indices[j]]->second.makeGoodForIterative();
				}
			}
		}
	}

	if ( cfg::get().expand_write_each_iteration ) {
		ofstream oftmp( getFilename(cfg::get().input_working_dir, Globals::iteration_no, "goodkmers", expand_iter_no ).data() );
		for ( hint_t n = 0; n < kmers.size(); ++n ) {
			if ( kmers[n]->second.isGoodForIterative() ) {
				oftmp << kmers[n]->first.str() << "\n>" << kmers[n]->first.start()
				      << "  cnt=" << kmers[n]->second.count << "  tql=" << (1-kmers[n]->second.totalQual) << "\n";
			}
		}
	}

	return res;
}

void HammerTools::PrintKMerResult( boost::iostreams::filtering_ostream & outf, const vector<KMerCount *> & kmers ) {
	for (vector<KMerCount *>::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
		outf << (*it)->first.start() << "\t"
			 << string(Globals::blob + (*it)->first.start(), K) << "\t"
			 << (*it)->second.count << "\t"
			 << (*it)->second.changeto << "\t"
			 << setw(8) << (*it)->second.totalQual << "\t";
		for (size_t i=0; i < K; ++i) outf << (*it)->second.qual[i] << " ";
		outf << "\n";
	}
}



bool HammerTools::CorrectOneRead( const vector<KMerCount*> & kmers, hint_t & changedReads,
		hint_t & changedNucleotides, hint_t readno, Read & r, size_t i ) {
	bool isGood = false;
	string seq (Globals::blob        + Globals::pr->at(readno).start(), Globals::pr->at(readno).size());
	string qual(Globals::blobquality + Globals::pr->at(readno).start(), Globals::pr->at(readno).size());
	const uint32_t read_size = Globals::pr->at(readno).size();
	PositionRead & pr = Globals::pr->at(readno);

	// create auxiliary structures for consensus
	vector<int> vA(read_size, 0), vC(read_size, 0), vG(read_size, 0), vT(read_size, 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);
	isGood = false;

	//cout << "  pr " << readno << " start=" << Globals::pr->at(readno).start() << "\tsize=" << Globals::pr->at(readno).size() << "\trevPos=" << Globals::revPos << "\trc=" << Globals::pr->at(readno).getRCBitString() << "\t" << Globals::pr->at(readno).getRCBit(0) << endl;
	//cout << string(Globals::blob + Globals::pr->at(readno).start(), Globals::pr->at(readno).size()) << endl;

	// getting the leftmost and rightmost positions of a solid kmer
	int left = read_size; int right = -1;
	bool changedRead = false;
	pair<int, hint_t> it = make_pair( -1, BLOBKMER_UNDEFINED );
	while ( (it = pr.nextKMerNo(it.first)).first > -1 ) {
		const PositionKMer & kmer = kmers[it.second]->first;
		const uint32_t pos = it.first;
		const KMerStat & stat = kmers[it.second]->second;
		//cout << "    pos=" << pos << "\trc=" << pr.getRCBit(pos) << "\tkmer=" << kmer.str() << endl;
		changedRead = changedRead || internalCorrectReadProcedure( r, pr, readno, seq, kmers, kmer, pos, stat, v, left, right, isGood, NULL, pr.getRCBit(pos) );
		//cout << "    internal ok " << " left=" << left << " right=" << right << endl;
	}

	int left_rev = 0; int right_rev = read_size-(int)K;

	if ( left <= right && left_rev <= right_rev ) {
		left = std::min(left, (int)read_size - left_rev - (int)K);
		right = std::max(right, (int)read_size - right_rev - (int)K);
	} else if ( left > right && left_rev <= right_rev ) {
		left = (int)read_size - left_rev - (int)K;
		right = (int)read_size - right_rev - (int)K;
	}

	// at this point the array v contains votes for consensus

	/*cout << "    consensus array:" << endl;
	for (size_t k=0; k<4; ++k) {
		for (size_t j=0; j<read_size; ++j) {
			cout << (char)((int)v[k][j] + (int)'0');
		}
		cout << endl;
	}*/

	size_t res = 0; // how many nucleotides have really changed?
	// find max consensus element
	for (size_t j=0; j<read_size; ++j) {
		char cmax = seq[j]; int nummax = 0;
		for (size_t k=0; k<4; ++k) {
			if (v[k][j] > nummax) {
				cmax = nucl(k); nummax = v[k][j];
			}
		}
		if (seq[j] != cmax) ++res;
		seq[j] = cmax;
	}

/*	cout << seq.data() << endl;
	cout << "  trimming to [" << left << ", " << right+K-1 << "]\tcurrent trim [" << r.ltrim() << ", " << r.rtrim() << "]" << endl;
	cout << "  before setsequence got read size=" << r.size() << ":\n    " << r.getSequenceString() << "\n    " << r.getPhredQualityString(64) << endl;
	cout << "  before setsequence got read size=" << r.size() << ":\n    " << r.getSequenceString() << "\n    " << r.getPhredQualityString(64) << endl;
	cout << "  before setsequence got read size=" << r.size() << ":\n    " << r.getSequenceString() << "\n    " << r.getPhredQualityString(64) << endl;*/
	r.setSequence(seq.data());
	r.trimLeftRight(left, right+(int)K-1);
	if ( left > 0 || right + K -1 < read_size ) changedRead = true;

	changedNucleotides += res;
	if (res > 0) ++changedReads;
	return isGood;
}


void HammerTools::CorrectReadFile( const string & readsFilename, const vector<KMerCount*> & kmers,
		hint_t & changedReads, hint_t & changedNucleotides, hint_t & readno, ofstream *outf_good, ofstream *outf_bad ) {
	ireadstream irs(readsFilename, cfg::get().input_qvoffset);
	VERIFY(irs.is_open());
	while (irs.is_open() && !irs.eof()) {
		Read r;
		irs >> r;
		size_t read_size = r.trimNsAndBadQuality(cfg::get().input_trim_quality);
		if (read_size < K) continue;
//		cout << "read:\t" << r.getSequenceString() << endl;
		if ( HammerTools::CorrectOneRead(kmers, changedReads, changedNucleotides, readno, r, 0) ) {
			r.print(*outf_good, cfg::get().input_qvoffset);
		} else {
			r.print(*outf_bad, cfg::get().input_qvoffset);
		}
		++readno;
	}
	irs.close();
}

void HammerTools::CorrectPairedReadFiles( const string & readsFilenameLeft, const string & readsFilenameRight,
		const vector<KMerCount*> & kmers, hint_t & changedReads, hint_t & changedNucleotides, hint_t readno_left_start, hint_t readno_right_start,
		ofstream * ofbadl, ofstream * ofcorl, ofstream * ofunpl, ofstream * ofbadr, ofstream * ofcorr, ofstream * ofunpr ) {
	ireadstream irsl(readsFilenameLeft, cfg::get().input_qvoffset);
	ireadstream irsr(readsFilenameRight, cfg::get().input_qvoffset);
	VERIFY(irsl.is_open());
	VERIFY(irsr.is_open());
	hint_t readno_left = readno_left_start;
	hint_t readno_right = readno_right_start;
	while (irsl.is_open() && !irsl.eof()) {
		Read l; Read r;
		irsl >> l; irsr >> r;
		size_t read_size_left = l.trimNsAndBadQuality(cfg::get().input_trim_quality);
		size_t read_size_right = r.trimNsAndBadQuality(cfg::get().input_trim_quality);
		if (read_size_left < K && read_size_right < K) {
			continue;
		}
		bool left_res = false;
		bool right_res = false;
		if (read_size_left >= K) {
			//cout << "  correcting lread " << readno_left << "\n" << l.getSequenceString() << "\n" << string(Globals::blob + Globals::pr->at(readno_left).start(), K) << endl;
			left_res = HammerTools::CorrectOneRead(kmers, changedReads, changedNucleotides, readno_left, l, 0);
			++readno_left;
		}
		if (read_size_right >= K) {
			//cout << "  correcting rread " << readno_right << "\n" << r.getSequenceString() << "\n" << string(Globals::blob + Globals::pr->at(readno_right).start(), K) << endl;
			right_res = HammerTools::CorrectOneRead(kmers, changedReads, changedNucleotides, readno_right, r, 0);
			++readno_right;
		}
		if (  !left_res ) l.print(*ofbadl, cfg::get().input_qvoffset);
		if ( !right_res ) r.print(*ofbadr, cfg::get().input_qvoffset);
		if (  left_res && !right_res ) l.print(*ofunpl, cfg::get().input_qvoffset);
		if ( !left_res &&  right_res ) r.print(*ofunpr, cfg::get().input_qvoffset);
		if (  left_res &&  right_res ) {
			l.print(*ofcorl, cfg::get().input_qvoffset);
			r.print(*ofcorr, cfg::get().input_qvoffset);
		}
	}
	irsl.close(); irsr.close();
}


hint_t HammerTools::CorrectAllReads() {
	// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
	hint_t changedReads = 0;
	hint_t changedNucleotides = 0;

	int correct_nthreads = min( cfg::get().correct_nthreads, cfg::get().general_max_nthreads );

	// TODO: make threaded read correction!
	TIMEDLN("Starting read correction in " << correct_nthreads << " threads (threaded correction does not work yet).");
	hint_t readno = 0;

	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		if (!cfg::get().input_paired) {
			ofstream ofgood(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "corrected.fastq").c_str());
			ofstream ofbad( HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "bad.fastq").c_str());
			HammerTools::CorrectReadFile(Globals::input_filenames[iFile], *Globals::kmers, changedReads, changedNucleotides, readno, &ofgood, &ofbad );
			TIMEDLN("  " << Globals::input_filenames[iFile].c_str() << " corrected.");
			// makes sense to change the input filenames for the next iteration immediately
			Globals::input_filenames[iFile] = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "corrected");
		} else {
			ofstream ofbadl(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "left.bad.fastq").c_str());
			ofstream ofcorr(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "right.corrected.fastq").c_str());
			ofstream ofbadr(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "right.bad.fastq").c_str());
			ofstream ofunpl(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "left.unpaired.fastq").c_str());
			ofstream ofunpr(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "right.unpaired.fastq").c_str());
			ofstream ofcorl(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "left.corrected.fastq").c_str());
			HammerTools::CorrectPairedReadFiles(Globals::input_filenames[iFile], Globals::input_filenames[iFile+1],
					*Globals::kmers, changedReads, changedNucleotides,
					Globals::input_file_blob_positions[iFile], Globals::input_file_blob_positions[iFile+1],
					&ofbadl, &ofcorl, &ofunpl, &ofbadr, &ofcorr, &ofunpr );
			TIMEDLN("  " << Globals::input_filenames[iFile].c_str() << " and " << Globals::input_filenames[iFile+1].c_str() << " corrected as a pair.");
			// makes sense to change the input filenames for the next iteration immediately
			Globals::input_filenames[iFile] = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "left.corrected.fastq");
			Globals::input_filenames[iFile+1] = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "reads", iFile, "right.corrected.fastq");
			++iFile;
		}
	}

	TIMEDLN("Correction done. Changed " << changedNucleotides << " bases in " << changedReads << " reads.");
	return changedReads;
}


void HammerTools::ReadKmerNosFromFile( const string & fname, vector<hint_t> *kmernos ) {
	kmernos->clear();
	boost::shared_ptr<FIStream> fis = FIStream::init(fname);
	hint_t pos;
	string buf;
	while (fis->fs.good()) {
		std::getline(fis->fs, buf);
		sscanf(buf.c_str(), "%lu", &pos);
		kmernos->push_back(pos);
	}
}

void HammerTools::ReadKmersAndNosFromFile( const string & fname, vector<KMerCount*> *kmers, vector<hint_t> *kmernos ) {
	kmernos->clear();
	kmers->clear();
	boost::shared_ptr<FIStream> fis = FIStream::init(fname);
	string buf;
	while (fis->fs.good()) {
		std::getline(fis->fs, buf);
		hint_t pos; int cnt; double qual;
		sscanf(buf.c_str(), "%lu\t%*s\t%u\t%lf", &pos, &cnt, &qual);
		kmernos->push_back(pos);
		kmers->push_back(new KMerCount(PositionKMer(pos), KMerStat(cnt, KMERSTAT_BAD, qual)));
	}
}

void HammerTools::ReadKmersWithChangeToFromFile( const string & fname, vector<KMerCount*> *kmers, vector<hint_t> *kmernos ) {
	kmernos->clear();
	kmers->clear();
	boost::shared_ptr<FIStream> fis = FIStream::init(fname);
	string buf;
	while (fis->fs.good()) {
		std::getline(fis->fs, buf);
		hint_t pos; int cnt; double qual; hint_t chg;
		sscanf(buf.c_str(), "%lu\t%*s\t%u\t%lu\t%lf", &pos, &cnt, &chg, &qual);
		kmernos->push_back(pos);
		kmers->push_back(new KMerCount(PositionKMer(pos), KMerStat(cnt, chg, qual)));
	}
}

string HammerTools::getFilename( const string & dirprefix, const string & suffix ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << suffix.data();
	return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, int iter_count, const string & suffix ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data();
	return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, const string & suffix, int suffix_num ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << suffix.data() << "." << suffix_num;
	return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data() << "." << suffix_num;
	return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num, const string & suffix2 ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data() << "." << suffix_num << "." << suffix2.data();
	return tmp.str();
}


bool HammerTools::internalCorrectReadProcedure( const Read & r, const PositionRead & pr, const hint_t readno, const string & seq, const vector<KMerCount*> & km,
		const PositionKMer & kmer, const uint32_t pos, const KMerStat & stat, vector< vector<int> > & v,
		int & left, int & right, bool & isGood, ofstream * ofs, bool revcomp ) {
	bool res = false;
	if (  stat.isGoodForIterative() || ( cfg::get().correct_use_threshold && stat.isGood() ) ) {
		isGood = true;
		//std::cout << "\t\t\tsolid " << kmer.str() << "\tpos=" << pos << endl;
		for (size_t j = 0; j < K; ++j) {
			if (!revcomp)
				v[dignucl(kmer[j])][pos + j]++;
			else
				v[complement(dignucl(kmer[j]))][K-1+pos-j]++;
		}
		if ((int) pos < left) left = pos;
		if ((int) pos > right) right = pos;
		//for (size_t j = 0; j < pos; ++j) std::cout << " ";
		//std::cout << kmer.str().data() << "\t" << pos << endl;

	} else {
		// if discard_only_singletons = true, we always use centers of clusters that do not coincide with the current center
		if (stat.change() && (cfg::get().bayes_discard_only_singletons
				|| km[stat.changeto]->second.isGoodForIterative()
				|| ( cfg::get().correct_use_threshold && stat.isGood() ))) {
			//std::cout << "\tchange to\n";
			//cout << "  kmer " << kmer.str() << " wants to change to " << km[stat.changeto]->first.str() << " pos=" << pos << endl;
			isGood = true;
			if ((int) pos < left) left = pos;
			if ((int) pos > right) right = pos;
			const PositionKMer & newkmer = km[stat.changeto]->first;

			for (size_t j = 0; j < K; ++j) {
				if (!revcomp)
					v[dignucl(newkmer[j])][pos + j]++;
				else
					v[complement(dignucl(newkmer[j]))][K-1+pos-j]++;
			}

			// pretty print the k-mer
			res = true;
			if (ofs != NULL) {
				for (size_t j = 0; j < pos; ++j)
					std::cout << " ";
				std::cout << newkmer.str().data();
			}
		}
	}
	return res;
}

void HammerTools::ProcessKmerHashFile( boost::iostreams::filtering_istream & inf, KMerNoHashMap & km ) {
	string buf; // a line contains two numbers, 1024 should be enough for everybody
	uint64_t pos; double prob;
	while (inf.good()) {
		std::getline(inf, buf);
		sscanf(buf.c_str(), "%lu\t%lf", &pos, &prob);
		KMerNo kmerno(pos, prob);
		KMerNoHashMap::iterator it_hash = km.find( kmerno );
		if ( it_hash == km.end() ) {
			KMerCount * kmc = new KMerCount( PositionKMer(kmerno.index), KMerStat(1, KMERSTAT_GOODITER, kmerno.errprob) );
			for (uint32_t j=0; j<K; ++j) {
				kmc->second.qual.set(j, Globals::blobquality[kmerno.index + j] - (char)cfg::get().input_qvoffset);
			}
			km.insert( make_pair( kmerno, kmc ) );
		} else {
			it_hash->second->second.count++;
			it_hash->second->second.totalQual *= kmerno.errprob;
			for (uint32_t j=0; j<K; ++j) {
				const int cur_qual = (int)it_hash->second->second.qual[j];
				if (cur_qual + (int)Globals::blobquality[kmerno.index + j] - (int)cfg::get().input_qvoffset < MAX_SHORT) {
					it_hash->second->second.qual.set(j, cur_qual + (short)Globals::blobquality[kmerno.index + j] - (short)cfg::get().input_qvoffset);
				} else {
					it_hash->second->second.qual.set(j, MAX_SHORT);
				}
			}
		}
	}
}

void HammerTools::PrintProcessedKmerHashFile( boost::iostreams::filtering_ostream & outf, hint_t & kmer_num, KMerNoHashMap & km ) {
	for (KMerNoHashMap::iterator it = km.begin(); it != km.end(); ++it) {
		outf << it->second->first.start() << "\t"
				<< string(Globals::blob + it->second->first.start(), K) << "\t"
				<< it->second->second.count << "\t"
				<< setw(8) << it->second->second.totalQual << "\t";
		for (size_t i=0; i < K; ++i) outf << it->second->second.qual[i] << " ";
		outf << "\n";
		delete it->second;
		++kmer_num;
	}
	km.clear();
}

FIStream::FIStream(const string & fname) : stdstream(fname, cfg::get().general_gzip ? (std::ios::in | std::ios::binary) : std::ios::in) {
	if (cfg::get().general_gzip) fs.push(boost::iostreams::gzip_decompressor());
	fs.push(stdstream);
}

FIStream::FIStream(const string & fname, bool input_output) : stdstream(fname,
		((!input_output && cfg::get().general_gzip) || (input_output && cfg::get().input_gzipped)) ? (std::ios::in | std::ios::binary) : std::ios::in) {
	if (cfg::get().general_gzip) fs.push(boost::iostreams::gzip_decompressor());
	fs.push(stdstream);
}

FIStream::~FIStream() { fs.reset(); }

boost::shared_ptr<FIStream> FIStream::init(const string & fname, bool input_output) {
	boost::shared_ptr<FIStream> p(new FIStream(fname, input_output));
	return p;
}

FOStream::FOStream(const string & fname) : stdstream(fname, cfg::get().general_gzip ? (std::ios::out | std::ios::binary) : std::ios::out) {
	if (cfg::get().general_gzip) fs.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
	fs.push(stdstream);
}

FOStream::FOStream(const string & fname, bool input_output) : stdstream(fname,
		((!input_output && cfg::get().general_gzip) || (input_output && cfg::get().input_gzipped)) ? (std::ios::out | std::ios::binary) : std::ios::out) {
	if (cfg::get().general_gzip) fs.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
	fs.push(stdstream);
}

FOStream::~FOStream() { fs.reset(); }

boost::shared_ptr<FOStream> FOStream::init(const string & fname, bool input_output) {
	boost::shared_ptr<FOStream> p(new FOStream(fname, input_output));
	return p;
}
