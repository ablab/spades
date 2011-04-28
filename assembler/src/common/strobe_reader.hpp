#ifndef STROBE_READER_HPP_
#define STROBE_READER_HPP_

using namespace std;

template<size_t cnt, typename T, typename TR>
class StrobeReader {
	vector<TR*> readers_;
public:
	typedef T ReadType;
	//	StrobeReader(const TR **readers) {
	//		for (size_t i = 0; i < cnt; ++i) {
	//			readers_[i] = readers[i];
	//		}
	//	}

	StrobeReader(const string filenames[]) {
		stringstream s;
		for (size_t i = 0; i < cnt; ++i) {
			readers_.push_back(new TR(filenames[i]));
		}
	}

	virtual ~StrobeReader() {
		close();
		for (size_t i = 0; i < cnt; ++i) {
			delete readers_[i];
		}
	}

	bool eof() const {
		for (size_t i = 0; i < cnt; ++i) {
			if (readers_[i]->eof()) {
				return true;
			}
		}
		return false;
	}

	StrobeReader& operator>>(vector<T>& v) {
		v.clear();
		T t;
		for (size_t i = 0; i < cnt; ++i) {
			(*readers_[i]) >> t;
			v.push_back(t);
		}
		return *this;
	}

	void reset() {
		for (size_t i = 0; i < cnt; ++i) {
			readers_[i]->reset();
		}
	}

	void close() {
		for (size_t i = 0; i < cnt; ++i) {
			readers_[i]->close();
		}
	}
};

template<typename T, typename TR>
struct MateReader {
	typedef StrobeReader<2, T, TR> type;
};

template<typename T, typename TR>
struct SingleReader {
	typedef StrobeReader<1, T, TR> type;
};

template<typename T, typename TR>
class CuttingReader {
	TR reader_;
	size_t cut_;
	size_t read_;

	CuttingReader(TR reader, size_t cut = -1) : reader_(reader), cut_(cut), read_(0) {}

	virtual ~CuttingReader() {
		close();
		delete reader_;
	}

	bool eof() const {
		return read_ == cut_ || reader_.eof();
	}

	CuttingReader& operator>>(T& v) {
		reader_ >> v;
		++read_;
		return *this;
	}

	void reset() {
		read_ = 0;
		reader_.reset();
	}

	void close() {
		reader_.close();
	}
};

#endif /* STROBE_READER_HPP_ */
