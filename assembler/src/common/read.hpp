/*
 * read.hpp
 *
 *  Created on: 29.03.2011
 *      Author: vyahhi
 */

#ifndef READ_HPP_
#define READ_HPP_

#include <string>
using namespace std;

class Read {
public:
	Read() {
		seq = qual = NULL;
		name = "undefined";
	}
private:
	Sequence *seq;
	Qual *qual;
	string name;
};

#endif /* READ_HPP_ */
