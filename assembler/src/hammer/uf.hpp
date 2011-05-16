/*
 * uf.hpp
 * Contains classes and procedures for working with .uf files
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
 
#ifndef UF_HPP
#define UF_HPP

#include<iostream>
#include<fstream>
#include<strstream>
#include<seq.hpp>

#define READ_BUFFER 256

using namespace std;

/**
 * A cluster of k-mers read from a uf file
 */
template<size_t ufsize_> class UFCluster {
  public:
	UFCluster() : center_(-1) { }
	
	void addSeq(const Seq<ufsize_> & seq) { vseq_.push_back(seq); }
	
	size_t center() const { return center_; }
	void setCenter(size_t center) { center_ = center; }
	
	const vector<Seq<ufsize_> > & vseq() const { return vseq_; }
	void setVSeq(const vector<Seq<ufsize_> > & vseq) { vseq_ = vseq; }
	
	const Seq<ufsize_> seq(size_t i) const { return vseq_[i]; }
	
	bool isValid() const { return (center_ >=0) && (vseq_.size() > 0) && (center_ < vseq_.size()); }
	
	void clear() { vseq_.clear(); center_ = -1; }
	size_t size() { return vseq_.size(); }
	
  private:
    vector<Seq<ufsize_> > vseq_;
    size_t center_;
};


/*
 * Read clusters from uf data (one by one)
 */
template<size_t size_> class iufstream {

public:
	iufstream(const string& filename) {
		filename_ = filename;
		is_open_ = open(filename);
	}

	virtual ~iufstream() {
		close();
	}

	bool is_open() const {
		return is_open_;
	}

	bool eof() const {
		return eof_;
	}

	static vector<UFCluster<size_> >* readAll(string filename, int cnt = -1) {
		iufstream irs(filename);
		assert(irs.is_open());
		vector<UFCluster<size_> >* res = new vector<UFCluster<size_> >();
		UFCluster<size_> ufc;
		while (cnt-- && irs.is_open() && !irs.eof()) {
			irs >> ufc;
			if (!ufc.isValid()) {
				cnt++;
				continue;
			}
			res->push_back(ufc);
		}
		irs.close();
		return res;
	}

	iufstream& operator>>(UFCluster<size_>  &ufc) {
		assert(is_open());
		assert(!eof());
		if (!is_open() || eof()) {
			return *this;
		}
		ufc.clear();
		ufc.setCenter(ufc_.center());
		ufc.setVSeq(ufc_.vseq());
		read_ahead(); // make actual read for the next result
		return *this;
	}

	void close() {
		if (is_open()) {
			ifs_.close();
			is_open_ = false;
		}
	}

	void reset() {
		ifs_.close();
		ifs_.open(filename_.c_str(), ifstream::in);
	}

private:
	std::string filename_;
	ifstream ifs_;
	kseq_t* seq_;
	bool is_open_;
	bool eof_;
	bool rtl_;
	UFCluster<size_> ufc_;
	
	/*
	 * open i's file with FASTQ reads,
	 * return true if it opened file, false otherwise
	 */
	bool open(string filename) {
		ifs_.open(filename.c_str(), ifstream::in); // STEP 2: open the file handler
		if (!ifs_.good()) {
			return false;
		}
		is_open_ = true;
		eof_ = false;
		read_ahead();
		return true;
	}

	void read_ahead() {
		assert(ifs_.good());
		
		char buffer[READ_BUFFER];
		istrstream ostr(buffer, READ_BUFFER);
		
		ufc_.clear();
		int index, count, total, dist;
		float freq, mult;
		char kmer[READ_BUFFER];
		char reason[READ_BUFFER];
		while (ifs_.good()) {
			ifs_.getline(buffer, READ_BUFFER);
			if (strlen(buffer) < size_) {// this is an empty line that may mark the end of a cluster
				if (ufc_.size() > 0) break;
				else continue;
			}
			sscanf(buffer, "%i\t%s\t%i\t%f\t%i\t%f\t%s\t%i", &index, kmer, &count, &freq, &total, &mult, reason, &dist);
			Seq<size_> seq(kmer);
			if (!strcmp(reason, "goodSingleton")) { // singleton; we drop them here
				if (ufc_.size() > 0) break;
				else continue;
			}
			ufc_.addSeq(seq);
			if (!strcmp(reason, "center")) ufc_.setCenter( ufc_.size() - 1 ); // center of a cluster
		}
		if (!ufc_.isValid()) {
			cout << index << " is invalid! size = " << ufc_.size() << "\n";
		}
		
		if (!ifs_.good()) eof_ = true;
	}
};

#endif 
