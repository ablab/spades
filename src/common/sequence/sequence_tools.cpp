#include "sequence_tools.hpp"

#include "utils/logger/logger.hpp"

#include "edlib/edlib.h"


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

