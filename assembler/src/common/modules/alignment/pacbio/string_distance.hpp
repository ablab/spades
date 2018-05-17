#include "edlib/edlib.h"

namespace pacbio {
static const int STRING_DIST_INF = 1e8;

inline int StringDistance(const string &a, const string &b) {
    int a_len = (int) a.length();
    int b_len = (int) b.length();
    int d = min(a_len / 3, b_len / 3);
    d = max(d, 10);
    if (a_len == 0 || b_len == 0) {
        d = max(a_len, b_len);
        if (d > 20) {
            INFO ("zero length path , lengths " << a_len << " and " << b_len);
            return STRING_DIST_INF;
        } else {
            return d;
        }
    }

    DEBUG(a_len << " " << b_len << " " << d);
    edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                                , {'R', 'A'}, {'R', 'G'}
                                                , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                                , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                                , {'M', 'A'}, {'M', 'C'}
                                                , {'S', 'C'}, {'S', 'G'}
                                                , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                                , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                                , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                                , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                                , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                                , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };
    edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len, b.c_str(), b_len
                                                   , edlib::edlibNewAlignConfig(2*d, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_DISTANCE,
                                                                         additionalEqualities, 36));
    int score = STRING_DIST_INF;
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        score = result.editDistance;
    }
    edlib::edlibFreeAlignResult(result);
    return score;
}

inline void SHWDistance(string &target, string &query, int max_score, vector<int> &positions, vector<int> &scores) {
    if (query.size() == 0){
        for (int i = 0; i < min(max_score, (int) target.size()); ++ i) {
            positions.push_back(i);
            scores.push_back(i + 1);
        }
        return;
    }
    if (target.size() == 0){
        if (query.size() <= max_score) {
		positions.push_back(0);
        	scores.push_back(query.size());
        }
        return;
    }
    VERIFY(target.size() > 0)
    edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                        , {'R', 'A'}, {'R', 'G'}
                                        , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                        , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                        , {'M', 'A'}, {'M', 'C'}
                                        , {'S', 'C'}, {'S', 'G'}
                                        , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                        , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                        , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                        , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                        , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                        , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };
    edlib::EdlibAlignResult result = edlib::edlibAlign(query.c_str(), (int) query.size(), target.c_str(), (int) target.size()
                                                   , edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_SHW_FULL, edlib::EDLIB_TASK_DISTANCE,
                                                                         additionalEqualities, 36));
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        positions.reserve(result.numLocations);
        scores.reserve(result.numLocations);
        for (int i = 0; i < result.numLocations; ++ i) {
            //INFO("Loc=" << result.endLocations[i] << " score=" << result.endScores[i]);
            if (result.endLocations[i] >= 0) {
                positions.push_back(result.endLocations[i]);
                scores.push_back(result.endScores[i]);
            }
        }
    }
    edlib::edlibFreeAlignResult(result);
}

inline int NWDistance(const string &a, const string &b, int max_score) {
    int a_len = (int) a.length();
    int b_len = (int) b.length();
    if (a_len == 0){
       if (b_len > max_score){
           return -1;
       } else{
           return b_len;
       }
    }
    if (b_len == 0){
       if (a_len > max_score){
           return -1;
       } else{
           return a_len;
       }
    }
    DEBUG(a_len << " " << b_len << " " << max_score);
    edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                                , {'R', 'A'}, {'R', 'G'}
                                                , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                                , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                                , {'M', 'A'}, {'M', 'C'}
                                                , {'S', 'C'}, {'S', 'G'}
                                                , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                                , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                                , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                                , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                                , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                                , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };

    edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len, b.c_str(), b_len
                                                   , edlib::edlibNewAlignConfig(max_score, edlib::EDLIB_MODE_NW, edlib::EDLIB_TASK_DISTANCE,
                                                                         additionalEqualities, 36));
    int score = -1;
    if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
        score = result.editDistance;
    }
    edlib::edlibFreeAlignResult(result);
    return score;
}

int SHWDistance2(const string &a, const string &b, int path_max_length_, int &end_pos) {
        //INFO("Before ed")
        int a_len = (int) a.length();
        int b_len = (int) b.length();
        VERIFY(a_len > 0);
        VERIFY(b_len > 0);
        edlib::EdlibEqualityPair additionalEqualities[36] = {{'U', 'T'}
                                                , {'R', 'A'}, {'R', 'G'}
                                                , {'Y', 'C'}, {'Y', 'T'}, {'Y', 'U'}
                                                , {'K', 'G'}, {'K', 'T'}, {'K', 'U'}
                                                , {'M', 'A'}, {'M', 'C'}
                                                , {'S', 'C'}, {'S', 'G'}
                                                , {'W', 'A'}, {'W', 'T'}, {'W', 'U'}
                                                , {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'B', 'U'}
                                                , {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'D', 'U'}
                                                , {'H', 'A'}, {'H', 'C'}, {'H', 'T'}, {'H', 'U'}
                                                , {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
                                                , {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}, {'N', 'U'} };
        edlib::EdlibAlignResult result = edlib::edlibAlign(a.c_str(), a_len, b.c_str(), b_len
                                                       , edlib::edlibNewAlignConfig(path_max_length_, edlib::EDLIB_MODE_SHW, edlib::EDLIB_TASK_DISTANCE,
                                                                             additionalEqualities, 36));
        int score = -1;
        if (result.status == edlib::EDLIB_STATUS_OK && result.editDistance >= 0) {
            if (result.numLocations > 0) {
                score = result.editDistance;
                end_pos = result.endLocations[0];
            } else {
                INFO("edlib: Strange")
            }
        }
        edlib::edlibFreeAlignResult(result);
        //INFO("After ed")
        //delete additionalEqualities;
        return score;
    }

}
