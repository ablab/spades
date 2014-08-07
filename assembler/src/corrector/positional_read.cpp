#include "positional_read.hpp"
using namespace std;

namespace corrector {

void position_description::update(const position_description &another) {
    for (size_t i = 0; i < MAX_VARIANTS; i++)
        votes_[i] += another.votes_[i];
    for (auto &ins : another.insertions_)
        insertions_[ins.first] += ins.second;
}

string position_description::str() const {
    stringstream ss;
    for (int i = 0; i < MAX_VARIANTS; i++) {
        ss << pos_to_var[i];
        ss << ": " << votes_[i] << "; ";
    }
    return ss.str();
}

size_t position_description::FoundOptimal(char current) const {
    size_t maxi = var_to_pos[(size_t) current];
    int maxx = votes_[maxi];
    for (size_t j = 0; j < MAX_VARIANTS; j++) {
        //1.5 because insertion goes _after_ match
        if (maxx < votes_[j] || (j == Variants::Insertion && maxx * 2 < votes_[j] * 3)) {
            maxx = votes_[j];
            maxi = j;
        }
    }
    return maxi;
}
void position_description::clear() {
    for (size_t i = 0; i < MAX_VARIANTS; i++) {
        votes_[i] = 0;
        insertions_.clear();
    }
}
}
;
