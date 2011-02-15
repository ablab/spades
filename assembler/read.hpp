#ifndef READ_HPP
#define READ_HPP

#include <string>

class Kmer {
	public:
		Kmer(const std::string &s);
		std::string to_string();
	private:
		long long val_; // max 32 bases
};

class SingleRead {
	public:
		SingleRead(const std::string &s);
		std::string to_string();
	private:
		long long val_[4]; // max 128 bases
};

class MateRead {
	public:
		MateRead(const SingleRead &s1, const SingleRead &s2);
		std::string to_string();
	private:
		SingleRead s1_;
		SingleRead s2_;
};

#endif
