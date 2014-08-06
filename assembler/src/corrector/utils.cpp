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
//Re: For this function wont fix. a) I want to save old NODE_ID, Cov and  UID;  b) while outputing this concrete contig we do not know number of processed in other threads, so we need to use one locking osequencestream for all threads..
//Changed string sum to  MakeContigId(int number, size_t length, double coverage, size_t id) {
// WTF: Extend osequencestream. Prepare patch and submit for review.

string ContigRenameWithLength(string name, size_t len) {
    vector < string > splitted = split(name, '_');
    if (splitted.size() >= 8 && splitted[0] == "NODE" && splitted[2] == "length") {
        string res = io::MakeContigId(stoi(splitted[1]), len, stod(splitted[5]), stoi(splitted[7]));
        return res;
    } else {
        return name;
    }
}

}
;
