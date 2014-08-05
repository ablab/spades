#include "utils.hpp"
#include "io/osequencestream.hpp"

#include <sstream>


using namespace std;

namespace corrector {
vector<string> split(const string &s, char delim) {
    vector < string > elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
// WTF: Get rid of this. Use osequencestream.
//Re: Wont fix. a) We to save old NODE and especially ID, b) while outputing this concrete contig we do not know number of processed in other threads, so we need to set NODE manually.
//Changed string sum to  MakeContigId(int number, size_t length, double coverage, size_t id) {

string ContigRenameWithLength(string name, size_t len) {
    vector < string > splitted = split(name, '_');
    if (splitted.size() == 8 && splitted[0] == "NODE" && splitted[2] == "length") {
        string res = io::MakeContigId(stoi(splitted[1]), len, stod(splitted[5]), stoi(splitted[7]));
        return res;
    } else {
        return name;
    }
}

// WTF: Why do you think that writing your own FASTA parser is a sane idea? Get rid of it, use what we're having in io
std::map<std::string, std::string> GetContigs(std::string filename) {
    std::map < std::string, std::string > res;
    std::ifstream istr(filename, std::ios_base::in);
    std::string seq = "";
    std::string cont_name = "";
    string line;
    stringstream ss;
    while (!istr.eof()) {
        line = "";
        istr >> line;
        if (line.length() > 0 && line[0] == '>') {
            if (cont_name != "") {

                res[cont_name] = ss.str();

            }
            ss.str().clear();
            ss.str("");
            ss.clear();
            cont_name = line.substr(1, line.length() - 1);
        } else {
            ss << line;
        }
    }
    if (cont_name != "") {
        res[cont_name] = ss.str();
    }
    istr.close();
    return res;
}

}
;
