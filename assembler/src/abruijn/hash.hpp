#ifndef HASH_HPP_
#define HASH_HPP_

template< typename T >
struct Hash {
public:
	unsigned int operator() (const T &seq) const;
};

template< typename T >
struct HashSym {
public:
	unsigned int operator() (const T &seq) const;
};

#endif /* HASH_HPP_ */
