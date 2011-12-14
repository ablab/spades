#ifndef IS_CORRUPTING_WRAPPER_HPP_
#define IS_CORRUPTING_WRAPPER_HPP_

namespace io {

class ISCorruptingWrapper: public DelegatingReaderWrapper<PairedRead> {
private:
	const size_t is_;
public:
	typedef PairedRead ReadType;

	explicit ISCorruptingWrapper(IReader<ReadType>& reader, size_t is) :
			DelegatingReaderWrapper(reader), is_(is) {
	}

	/* virtual */
	ISCorruptingWrapper& operator>>(ReadType& read) {
		(*reader_) >> read;
		read = PairedRead(read.first(), read.second(), is_);
		return *this;
	}
};

}

#endif /* IS_CORRUPTING_WRAPPER_HPP_ */
