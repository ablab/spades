#include "QcException.hpp"

QcException::QcException(const std::string m) {
	this->message = m;
};

QcException::~QcException() throw (){};

const char* QcException::what() const throw () {
	return message.c_str();
};
