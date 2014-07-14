#include "utils.hpp"
#include <sstream>

using namespace std;

namespace corrector {
vector<string> split(const string &s, char delim) {
	vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

string ContigRenameWithLength( string name, size_t len) {
	vector<string> splitted = split(name, '_');
	if (splitted.size() == 6 && splitted[0] == "NODE" && splitted[2] == "length") {
		splitted[3] = std::to_string(len);
		string res = "";
		for (size_t i = 0; i < splitted.size() - 1; i++)
			res += splitted[i] + "_";
		res += splitted[splitted.size() - 1];
		return res;
	} else {
		return name;
	}
}


};
