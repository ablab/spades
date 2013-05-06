#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <exception>
#include <string>

class QcException : public std::exception {
public:
	QcException(const std::string m);
	~QcException() throw ();

	const char* what() const throw ();

private:
	std::string message;
};

#endif /* EXCEPTION_H_ */
