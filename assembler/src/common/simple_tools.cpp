
#include "simple_tools.hpp"
#include "nucl.hpp"

string Reverse(const string &s) {
	return string(s.rbegin(), s.rend());
}

string Complement(const string &s) {
	string res = s;
	for (size_t i = 0; i < s.size(); i++) {
		if (res[i] != 'N' && res[i] != 'n') {
			res[i] = nucl_complement(res[i]);
		}
	}
	return res;
}

string ReverseComplement(const string &s) {
	return Complement(Reverse(s));
}
