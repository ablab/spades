#include "read.hpp"
#include "include.hpp"
using namespace std;

// FIXME: EVERYWHERE: USE SPACES, NOT TABS! FIX ALL THE CODING STYLE PROBLEMS EVERYWHERE

namespace corrector {

void position_description::update(const position_description &another) {
    for (size_t i = 0; i < MAX_VARIANTS; i++)
        votes[i] += another.votes[i];
    for (auto &ins : another.insertions)
        insertions[ins.first] += ins.second;
}
string position_description::str() const {
    stringstream ss;
    for (int i = 0; i < MAX_VARIANTS; i++) {
        ss << pos_to_var[i];
        ss << ": " << votes[i] << "; ";

    }
    return ss.str();
}

size_t position_description::FoundOptimal(char current) const {
    size_t maxi = var_to_pos[(size_t) current];
    int maxx = votes[maxi];
    for (size_t j = 0; j < MAX_VARIANTS; j++) {
        //1.5 because insertion goes _after_ match
        //std::min_element
        if (maxx < votes[j] || (j == Variants::Insertion && maxx * 2 < votes[j] * 3)) {
            maxx = votes[j];
            maxi = j;
        }
    }
    return maxi;
}
void position_description::clear() {
    for (size_t i = 0; i < MAX_VARIANTS; i++) {
        votes[i] = 0;
        insertions.clear();
    }
}
}
;
