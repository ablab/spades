#include "read.hpp"
#include "include.hpp"
#include "sam_reader.hpp"
using namespace std;

// FIXME: EVERYWHERE: USE SPACES, NOT TABS! FIX ALL THE CODING STYLE PROBLEMS EVERYWHERE

namespace corrector {
    // WTF: const
	bool MappedSamStream::eof() {
		return eof_;
	}
    // WTF: const
    bool MappedSamStream::is_open() {
    	return is_open_;
    }

    MappedSamStream& MappedSamStream::operator>>(SingleSamRead& read) {
	  if (!is_open_ || eof_)
		  return *this;
	  read.set_data(seq_);
      // WTF: How samread can throw if this is a C library?
	  try {
		  int tmp = samread(reader_, seq_);
		  eof_ = (0 >= tmp);
	  } catch (std::exception e) {
		  WARN("Error in sam file" << filename_ << "skipping the rest");
          // WTF: 1. Why you're resetting eof_ here?
          // WTF: 2. Why you're using 0 for bool variable? Forgot about true and false?
		  eof_ = 0;
	  }


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
        // WTF: Make sure the message here is readable
    	VERIFY_MSG (r1.GetName() == r2.GetName(), r1.GetName() + " " + r2.GetName());
    	read.pair(r1,r2);
        return *this;
    }
    // WTF: this is getter. Name it properly. Why doesn't it const?
    bam_header_t* MappedSamStream::ReadHeader(){
    	return reader_->header;
    }

    // WTF: const
    string MappedSamStream::get_contig_name(int i){
    	VERIFY(i < reader_->header->n_targets);
    	return (reader_->header->target_name[i]);
    }
    // WTF: What's about seq_? Leaking?
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
            // WTF: Have you completely forgotten about the coding style?
		{
            // WTF: Make you messages readable
		   WARN( "Fail to open SAM/BAM file " << filename_ << " ,skipping");
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
