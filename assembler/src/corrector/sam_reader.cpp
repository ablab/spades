#include "read.hpp"
#include "include.hpp"
#include "sam_reader.hpp"
using namespace std;

namespace corrector {
	bool MappedSamStream::eof() {
		return eof_;
	}
    bool MappedSamStream::is_open() {
    	return is_open_;
    }

    MappedSamStream& MappedSamStream::operator>>(SingleSamRead& read) {
	  if (!is_open_ || eof_)
		  return *this;
	  read.set_data(seq_);
	  int tmp = samread(reader_, seq_);
	  eof_ = (0 >= tmp);
	  return *this;
    }

    MappedSamStream& MappedSamStream::operator >> (PairedSamRead& read){
    	TRACE("starting process paired read");
    	SingleSamRead r1;
    	MappedSamStream::operator >> (r1);
    	SingleSamRead r2;
    	MappedSamStream::operator >> (r2);
    	TRACE(r1.GetSeq());
    	TRACE(r2.GetSeq());
    	TRACE(r1.GetName());
    	VERIFY_MSG (r1.GetName() == r2.GetName(), r1.GetName() + " " + r2.GetName());
    	read.pair(r1,r2);
        return *this;
    }
    bam_header_t* MappedSamStream::ReadHeader(){
    	return reader_->header;
    }

    string MappedSamStream::get_contig_name(int i){
    	VERIFY(i < reader_->header->n_targets);
    	return (reader_->header->target_name[i]);
    }
    void MappedSamStream::close() {
    	samclose(reader_);
    	is_open_ = false;
    	eof_ = true;
    }

    void MappedSamStream::reset() {
    	close();
    	open();
    }

    void MappedSamStream::open() {
		if ((reader_ = samopen(filename_.c_str(), "r", NULL)) == NULL)
		{
		   cerr << "Fail to open SAM/BAM file " << filename_ << endl;
		   is_open_ = false;
		   eof_ = true;
		} else {
			is_open_ = true;
			int tmp = samread(reader_, seq_);
			eof_ = (0 >= tmp);
		}
	}
    io::ReadStreamStat MappedSamStream::get_stat() const {
    	return io::ReadStreamStat();
    }


}
