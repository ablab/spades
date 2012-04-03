#pragma once

#include "io/reader.hpp"
#include "io/delegating_reader_wrapper.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/filtering_reader_wrapper.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"

namespace io {
//todo refactor, and maybe merge them once again
class EasyReader : public DelegatingReaderWrapper<io::SingleRead> {
	explicit EasyReader(const EasyReader& reader);
	void operator=(const EasyReader& reader);

	Reader raw_reader_;
//	FilteringReaderWrapper<ReadType> filtered_reader_;
	CarefulFilteringReaderWrapper<io::SingleRead> filtered_reader_;
	RCReaderWrapper<io::SingleRead> rc_reader_;

public:
  explicit EasyReader(const io::SingleRead::FilenameType& filename,
                  OffsetType offset_type = PhredOffset)
      : raw_reader_(filename, offset_type),
        filtered_reader_(raw_reader_),
        rc_reader_(filtered_reader_) {
	  Init(rc_reader_);
  }

  /*
   * Default destructor.
   */
  /* virtual */ ~EasyReader() {
  }

};

class PairedEasyReader
	: public DelegatingReaderWrapper<io::PairedRead>
{
	scoped_ptr<IReader<io::PairedRead>> raw_reader_;
	CarefulFilteringReaderWrapper<io::PairedRead> filtered_reader_;
	RCReaderWrapper<io::PairedRead> rc_reader_;

public:
  PairedEasyReader(const io::PairedRead::FilenamesType& filenames,
                  size_t insert_size,
                  bool change_read_order = false,
                  OffsetType offset_type = PhredOffset)
      : raw_reader_(new SeparateReader(filenames, insert_size, change_read_order, offset_type))
  	  , filtered_reader_(*raw_reader_)
  	  , rc_reader_(filtered_reader_)
  {
	  Init(rc_reader_);
  }

  PairedEasyReader(const std::string& filename,
                  size_t insert_size,
                  bool change_read_order = false,
                  OffsetType offset_type = PhredOffset)
      : raw_reader_(new MixedReader(filename, insert_size, change_read_order, offset_type))
  	  , filtered_reader_(*raw_reader_)
  	  , rc_reader_(filtered_reader_)
  {
	  Init(rc_reader_);
  }
};

class PairedPureEasyReader : public DelegatingReaderWrapper<io::PairedRead> {
//	explicit PairedPureEasyReader(const PairedPureEasyReader& reader);
//	void operator=(const PairedPureEasyReader& reader);

	scoped_ptr<IReader<io::PairedRead>>  raw_reader_;
	CarefulFilteringReaderWrapper<io::PairedRead> filtered_reader_;

public:
  explicit PairedPureEasyReader(const io::PairedRead::FilenamesType& filenames,
                  size_t insert_size,
                  bool change_read_order = false,
                  OffsetType offset_type = PhredOffset)
      : raw_reader_(new SeparateReader(filenames, insert_size, change_read_order, offset_type, false)),
        filtered_reader_(*raw_reader_){
	  Init(filtered_reader_);
  }

  explicit PairedPureEasyReader(const std::string& filename,
                  size_t insert_size,
                  bool change_read_order = false,
                  OffsetType offset_type = PhredOffset)
      : raw_reader_(new MixedReader(filename, insert_size, change_read_order, offset_type, false)),
        filtered_reader_(*raw_reader_){
	  Init(filtered_reader_);
  }
  /*
   * Default destructor.
   */
  /* virtual */ ~PairedPureEasyReader() {
  }

};

class PureEasyReader : public DelegatingReaderWrapper<io::SingleRead> {
	explicit PureEasyReader(const PureEasyReader& reader);
	void operator=(const PureEasyReader& reader);

	Reader raw_reader_;
	CarefulFilteringReaderWrapper<io::SingleRead> filtered_reader_;

public:
  explicit PureEasyReader(const io::SingleRead::FilenameType& filename,
                  OffsetType offset_type = PhredOffset)
      : raw_reader_(filename, offset_type),
        filtered_reader_(raw_reader_){
	  Init(filtered_reader_);
  }

  /*
   * Default destructor.
   */
  /* virtual */ ~PureEasyReader() {
  }
};
}
