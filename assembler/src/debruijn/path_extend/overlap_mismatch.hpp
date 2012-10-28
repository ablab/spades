/*
 * overlap_mismatch.hpp
 *
 *  Created on: Oct 26, 2012
 *      Author: andrey
 */

#ifndef OVERLAP_MISMATCH_HPP_
#define OVERLAP_MISMATCH_HPP_

namespace path_extend {

using namespace debruijn_graph;
using namespace std;

class ContigOverlapMismatchCorrecter {

private:

    struct NuclMismatchInfo {

        size_t position;

        char new_symbol;
    };

    struct OverlapMismatchInfo {

        vector<NuclMismatchInfo> overlap_mismatches;

    };

    conj_graph_pack& gp_;

    map <size_t, vector <OverlapMismatchInfo> > nucl_map_;

    vector<size_t> DiffPos(const Sequence& s1, const Sequence& s2) const {
        VERIFY(s1.size() == s2.size());
        vector < size_t > answer;
        for (size_t i = 0; i < s1.size(); ++i)
            if (s1[i] != s2[i])
                answer.push_back(i);
        return answer;
    }

    size_t HammingDistance(const Sequence& s1, const Sequence& s2) const {
        VERIFY(s1.size() == s2.size());
        size_t dist = 0;
        for (size_t i = 0; i < s1.size(); ++i)
            if (s1[i] != s2[i])
                dist++;
        return dist;
    }

    OverlapMismatchInfo ProcessOverlap(EdgeId left, EdgeId right, int gap) {

    }

    vector<OverlapMismatchInfo> ProcessPath(const BidirectionalPath& path) {
        vector <OverlapMismatchInfo> overlaps(path.Size());

        for (size_t i = 1; i < path.Size(); ++i) {

            if (path.GapAt(i) < (int) gp_.g.k()) {
                int l = path.GapAt(i);
                EdgeId left = path.At(i - 1);
                EdgeId right = path.At(i);

                if (HammingDistance(gp_.g.EdgeNucls(left).Subseq(gp_.g.length(left) + gp_.g.k() - l), gp_.g.EdgeNucls(right).Subseq(0, l)) > 0) {

                }
            }
        }

        return overlaps;
    }


public:

};

}



#endif /* OVERLAP_MISMATCH_HPP_ */
