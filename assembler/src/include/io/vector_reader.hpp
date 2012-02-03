#pragma once
namespace io {

/**
 * Use vector<T> as input-stream with operator>>(T& t)
 */
template <typename T>
class VectorReader : public io::IReader<T> {
       std::vector<T> data_;
       size_t pos_;
       bool closed_;
public:
       VectorReader(const std::vector<T>& data) : data_(data), pos_(0), closed_(false) {

       }

       VectorReader(const T& item) : data_({item}), pos_(0), closed_(false) {

       }

       virtual ~VectorReader() {

       }

       virtual bool eof() /*const */{
               return pos_ == data_.size();
       }

       VectorReader<T>& operator>>(T& t) {
               t = data_[pos_++];
               return *this;
       }

       void close() {
               closed_ = true;
       }

       virtual bool is_open() /*const */{
               return !closed_;
       }

       void reset() {
               pos_ = 0;
       }

};

}
