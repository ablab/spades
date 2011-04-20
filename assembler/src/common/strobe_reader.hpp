#ifndef STROBE_READER_HPP_
#define STROBE_READER_HPP_

using namespace std;

template<size_t cnt, typename T, typename TR>
class StrobeReader {
	vector<TR*> readers_;
public:
	//	StrobeReader(const TR **readers) {
	//		for (size_t i = 0; i < cnt; ++i) {
	//			readers_[i] = readers[i];
	//		}
	//	}

	StrobeReader(const string *filenames) {
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

template<size_t count, typename T, typename TR>
class SimpleReaderWrapper {
private:
	StrobeReader<count, T, TR> &inner_reader_;
	vector<T> result_;
	int current_;
public:
	SimpleReaderWrapper(StrobeReader<count, T, TR> &reader) :
		inner_reader_(reader), current_(0) {
	}

	bool eof() const {
		return inner_reader_.eof();
	}

	SimpleReaderWrapper& operator>>(T& v) {
		if (current_ == 0) {
			inner_reader_ >> result_;
		}
		v = result_[current_];
		current_++;
		if (current_ == count)
			current_ = 0;
		return *this;
	}

	void reset() {
		current_ = 0;
		inner_reader_.reset();
	}

	void close() {
		inner_reader_.close();
	}
};

#endif /* STROBE_READER_HPP_ */
