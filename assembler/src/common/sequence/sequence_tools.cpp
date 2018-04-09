#include "edlib/edlib.h"
#include "sequence_tools.hpp"

int StringDistance(const std::string &a, const std::string &b) {
    int a_len = (int) a.length();
    int b_len = (int) b.length();
    int d = std::min(a_len / 3, b_len / 3);
    d = std::max(d, 10);
    if (a_len == 0 || b_len == 0) {
        if (d > 10) {
            // TRACE("zero length path , lengths " << a_len << " and " << b_len);
            return std::numeric_limits<int>::max();
        } else {
            return d;
        }
    }

    // DEBUG(a_len << " " << b_len << " " << d);
    edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len,
                                                       b.c_str(), b_len,
                                                       edlib::edlibNewAlignConfig(2*d, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_DISTANCE,
                                                                                  NULL, 0));
    int score = std::numeric_limits<int>::max();
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        score = result.editDistance;
    }
    edlib::edlibFreeAlignResult(result);
    return score;
}
