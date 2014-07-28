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


std::map<std::string, std::string> GetContigs(std::string filename) {
	std::map<std::string, std::string> res;
	std::ifstream istr(filename, std::ios_base::in);
    std::string seq = "";
    std::string cont_name = "";
    string line;
    stringstream ss;
    while(! istr.eof()) {
    	line = "";
    	istr >> line;
        if (line.length() > 0 && line[0] == '>') {
            if (cont_name != "") {

                res[cont_name] = ss.str();

            }
            ss.clear();
            cont_name = line.substr(1, line.length() - 1);
        } else {
            ss << line;
        }
    }
    if (cont_name != "") {
        res[cont_name] = ss.str();

    }
    return res;
}
void PutContig(std::string full_path, std::string contig_name, std::string contig_seq) {
	std::ofstream oss(full_path, std::ios_base::out);
	oss << ">" << contig_name << std::endl;
	size_t cur = 0;
	while (cur < contig_seq.size()) {
		oss << contig_seq.substr(cur, 60) << std::endl;
	    cur += 60;
	}
}


};
