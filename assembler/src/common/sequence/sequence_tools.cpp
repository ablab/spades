#include "sequence_tools.hpp"

#include "utils/logger/logger.hpp"

#include "edlib/edlib.h"

#include <map>

static std::map<std::string, std::string> protein_canonical = {
   {"A","GCT"},{"R","CGT"},{"N","AAT"},{"D","GAT"},{"C","TGT"},{"Q","CAA"},
   {"E","GAA"},{"G","GGT"},{"H","CAT"},{"I","ATT"},{"M","ATG"},{"L","TTA"},
   {"K","AAA"},{"F","TTT"},{"P","CCT"},{"S","TCT"},{"T","ACT"},{"W","TGG"},
   {"Y","TAT"},{"V","GTT"},{"*","TAG"}
};

static std::map<std::string, std::string> nuc_canonical = {
   {"GCT","A"}, {"GCC","A"}, {"GCA","A"}, {"GCG","A"}, {"CGT","R"}, {"CGC","R"}, {"CGA","R"}, {"CGG","R"}, 
   {"AGA","R"}, {"AGG","R"}, {"AAT","N"}, {"AAC","N"}, {"GAT","D"}, {"GAC","D"}, {"TGT","C"}, {"TGC","C"}, 
   {"CAA","Q"}, {"CAG","Q"}, {"GAA","E"}, {"GAG","E"}, {"GGT","G"}, {"GGA","G"}, {"GGC","G"}, {"GGG","G"}, 
   {"CAT","H"}, {"CAC","H"}, {"ATT","I"}, {"ATA","I"}, {"ATC","I"}, {"TTA","L"}, {"TTG","L"}, {"CTT","L"}, 
   {"CTC","L"}, {"CTA","L"}, {"CTG","L"}, {"AAA","K"}, {"AAG","K"}, {"ATG","M"}, {"TTT","F"}, {"TTC","F"}, 
   {"CCT","P"}, {"CCA","P"}, {"CCG","P"}, {"CCC","P"}, {"TCT","S"}, {"TCC","S"}, {"TCA","S"}, {"TCG","S"}, 
   {"AGT","S"}, {"AGC","S"}, {"ACT","T"}, {"ACC","T"}, {"ACA","T"}, {"ACG","T"}, {"TGG","W"}, {"TAT","Y"}, 
   {"TAC","Y"}, {"GTT","V"}, {"GTA","V"}, {"GTG","V"}, {"GTC","V"}, {"TAA","*"}, {"TGA","*"}, {"TAG","*"}
};

static std::map<char, int> proteins = {
        {'A', 0}, {'R', 1}, {'N', 2},
        {'D', 3}, {'C', 4}, {'Q', 5},
        {'E', 6}, {'G', 7}, {'H', 8},
        {'I', 9}, {'L', 11}, {'K', 12},
        {'M', 12}, {'F', 13}, {'P', 14},
        {'S', 15}, {'T', 16}, {'W', 17},
        {'Y', 18}, {'V', 19}, {'B', 20},
        {'Z', 21}, {'X', 22}, {'*', 23}
};

std::string ConvertProtein2CanonicalNuc(const std::string &seq){
    std::string res;
    for (auto c: seq) {
        if (protein_canonical.count(std::string(1, c)) > 0) {
            res += protein_canonical[std::string(1, c)];
        } else {
            res += "ATG";
        }
    }
    res += protein_canonical["*"];
    return res;
}

std::string ConvertNuc2CanonicalNuc(const std::string &seq){
    std::string res;
    for (size_t i = 0; i + 2 < seq.size(); i += 3) {
        std::string triplet = seq.substr(i, 3);
        if (nuc_canonical.count(triplet) > 0) {
            res += protein_canonical[nuc_canonical[triplet]];
        } else {
            res += "ATG";
        }
    }
    return res;
}

std::string ConvertNuc2Protein(const std::string &seq){
    std::string res;
    for (size_t i = 0; i + 2 < seq.size(); i += 3) {
        std::string triplet = seq.substr(i, 3);
        if (nuc_canonical.count(triplet) > 0) {
            res += nuc_canonical[triplet];
        } else {
            res += "M";
        }
    }
    return res;
}

bool EqualTriplets(const std::string &a, const std::string &b) {
    if (a.size() < 3 or b.size() < 3){
        return false;
    }
    return nuc_canonical[a] == nuc_canonical[b];
}

int BinaryScoreTriplets(const std::string &a, const std::string &b) {
    if (a.size() < 3 or b.size() < 3){
        return 1;
    }
    return nuc_canonical[a] == nuc_canonical[b] ? 0 : 1;
}

int ProteinStringDistance(const std::string &a, const std::string &b, int max_score) {
    int score = std::numeric_limits<int>::max();
    if (a.size() % 3 != 0 || b.size() % 3 != 0) return score;

    std::string p_a = ConvertNuc2Protein(a);
    std::string p_b = ConvertNuc2Protein(b);
    int a_len = (int) p_a.length();
    int b_len = (int) p_b.length();
    int d = std::min(a_len / 3, b_len / 3);
    d = std::max(d, 10);
    if (a_len == 0 || b_len == 0) {
        if (d > 10) {
            DEBUG("zero length path , lengths " << a_len << " and " << b_len);
            return std::numeric_limits<int>::max();
        } else {
            return std::max(a_len, b_len);
        }
    }
    if (max_score == -1) {
        max_score = 2 * d;
    }
    edlib::EdlibAlignResult result = edlib::edlibAlign(p_a.c_str(), a_len,
                                     p_b.c_str(), b_len,
                                     edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_DISTANCE, NULL, 0));
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        score = result.editDistance;
    }
    edlib::edlibFreeAlignResult(result);
    return score;
}


int StringDistance(const std::string &a, const std::string &b, int max_score) {
    int a_len = (int) a.length();
    int b_len = (int) b.length();
    int d = std::min(a_len / 3, b_len / 3);
    d = std::max(d, 10);
    if (a_len == 0 || b_len == 0) {
        if (d > 10) {
            DEBUG("zero length path , lengths " << a_len << " and " << b_len);
            return std::numeric_limits<int>::max();
        } else {
            return std::max(a_len, b_len);
        }
    }
    if (max_score == -1) {
        max_score = 2 * d;
    }
    edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len,
                                     b.c_str(), b_len,
                                     edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_DISTANCE, NULL, 0));
    int score = std::numeric_limits<int>::max();
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        score = result.editDistance;
    }
    edlib::edlibFreeAlignResult(result);
    return score;
}


void SHWDistanceExtended(const std::string &target, const std::string &query, int max_score, std::vector<int> &positions, std::vector<int> &scores) {
    if (query.size() == 0) {
        for (int i = 0; i < std::min(max_score, (int) target.size()); ++ i) {
            positions.push_back(i);
            scores.push_back(i + 1);
        }
        return;
    }
    if (target.size() == 0) {
        if ((int) query.size() <= max_score) {
            positions.push_back(0);
            scores.push_back((int) query.size());
        }
        return;
    }
    VERIFY(target.size() > 0)
    edlib::EdlibAlignResult result = edlib::edlibAlign(query.c_str(), (int) query.size(), target.c_str(), (int) target.size()
                                     , edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_SHW_EXTENDED, edlib::EDLIB_TASK_DISTANCE, NULL, 0));
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        positions.reserve(result.numLocations);
        scores.reserve(result.numLocations);
        for (int i = 0; i < result.numLocations; ++ i) {
            if (result.endLocations[i] >= 0) {
                positions.push_back(result.endLocations[i]);
                scores.push_back(result.endScores[i]);
            }
        }
    }
    edlib::edlibFreeAlignResult(result);
}

int SHWDistance(const std::string &a, const std::string &b, int max_score, int &end_pos) {
    int a_len = (int) a.length();
    int b_len = (int) b.length();
    VERIFY(a_len > 0);
    VERIFY(b_len > 0);
    edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len, b.c_str(), b_len
                                     , edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_SHW, edlib::EDLIB_TASK_DISTANCE, NULL, 0));
    int score = std::numeric_limits<int>::max();
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        if (result.numLocations > 0) {
            score = result.editDistance;
            end_pos = result.endLocations[0];
        }
    }
    edlib::edlibFreeAlignResult(result);
    return score;
}

