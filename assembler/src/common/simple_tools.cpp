#include "simple_tools.hpp"
#include "nucl.hpp"
string Reverse(const string &s) {
	size_t length = s.length();
	string result(length, 0);
	for (size_t i = 0; i < length; i++) {
		result[length - 1 - i] = s[i];
	}
	return result;
}

string Complement(const string &s) {
	size_t length = s.length();
	string result(length, 0);
	for (size_t i = 0; i < length; i++) {
		if (s[i] == 'N')
			result[i] = s[i];
		else
			result[i] = nucl_complement(s[i]);
	}
	return result;
}

string ReverseComplement(const string &s) {
	size_t length = s.length();
	string result(length, 0);
	for (size_t i = 0; i < length; i++) {
		if (s[i] == 'N')
			result[length - 1 - i] = s[i];
		else
			result[length - 1 - i] = nucl_complement(s[i]);
	}
	return result;
}
