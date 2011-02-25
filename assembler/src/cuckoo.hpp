/*
 * cuckoo.hpp
 *
 *  Created on: 25.02.2011
 *      Author: vyahhi
 */

#ifndef CUCKOO_HPP_
#define CUCKOO_HPP_

template <class key, class value, class hash = hash<key> >
class cuckoo {
public:
	cuckoo() {};
	virtual ~cuckoo() {};
};

#endif /* CUCKOO_HPP_ */
