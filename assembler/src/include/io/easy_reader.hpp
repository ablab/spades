#pragma once

#include "io/reader.hpp"
#include "io/delegating_reader_wrapper.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/filtering_reader_wrapper.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"

namespace io {

template<typename ReadType>
class EasyReader : public DelegatingReaderWrapper<ReadType> {
	explicit EasyReader(const EasyReader<ReadType>& reader);
	void operator=(const EasyReader<ReadType>& reader);

	Reader<ReadType> raw_reader_;
//	FilteringReaderWrapper<ReadType> filtered_reader;
	CarefulFilteringReaderWrapper<ReadType> filtered_reader_;
	RCReaderWrapper<ReadType> rc_reader_;

public:
  explicit EasyReader(const typename ReadType::FilenameType& filename,
                  size_t distance = 0,
                  OffsetType offset_type = PhredOffset)
      : raw_reader_(filename, distance, offset_type),
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

}
