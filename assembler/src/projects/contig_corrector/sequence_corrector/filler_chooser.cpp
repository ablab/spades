#include "filler_chooser.hpp"
#include "common/assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include "utils/logger/logger.hpp"

#include "edlib/edlib.h"

#include <utility>
#include <algorithm>
#include <sstream>

namespace sequence_corrector {

namespace {

struct AlignResultGuard : edlib::EdlibAlignResult {
    using Base = edlib::EdlibAlignResult;

    AlignResultGuard(Base res)
        : Base(res)
    {}

    AlignResultGuard(const AlignResultGuard & ) = delete;
    AlignResultGuard & operator=(const AlignResultGuard &) = delete;
    AlignResultGuard(AlignResultGuard &&) = delete;
    AlignResultGuard & operator=(AlignResultGuard &&) = delete;

    ~AlignResultGuard() {
        edlib::edlibFreeAlignResult(*this);
    }
};

edlib::EdlibAlignResult EditDistance(const std::string &a, const std::string &b, edlib::EdlibAlignTask mode = edlib::EDLIB_TASK_DISTANCE) {
    auto d = std::min(a.size() / 3, b.size() / 3);
    d = 2 * std::max(d, 10ul);
    // int d = -1;

    auto config = edlib::edlibNewAlignConfig(static_cast<int>(d), edlib::EDLIB_MODE_NW, mode, nullptr, 0);
    return edlib::edlibAlign(a.data(), a.size(), b.data(), b.size(), config);
}

struct AlignmentIterator {
    const unsigned char * alignment;
    int alignment_length;
    int align_pos;
    int s1_pos;
    int s2_pos;
    AlignmentIterator(AlignResultGuard const & alignment)
        : alignment(alignment.alignment)
        , alignment_length(alignment.alignmentLength)
        , align_pos(0)
        , s1_pos(0)
        , s2_pos(0)
    {}
    void operator ++ () noexcept {
        VERIFY(!IsEnd());
        if (alignment[align_pos] == EDLIB_EDOP_MATCH || alignment[align_pos] == EDLIB_EDOP_MISMATCH) {
            ++s1_pos;
            ++s2_pos;
        } else {
            if (alignment[align_pos] == EDLIB_EDOP_INSERT) 
                ++s1_pos;
            else
                ++s2_pos;
        }
        ++align_pos;
    }

    bool IsEnd() const noexcept {
        return alignment_length <= align_pos;
    }

    unsigned char Operation() const noexcept {
        return alignment[align_pos];
    }

    bool IsMatch() const noexcept {
        return Operation() == EDLIB_EDOP_MATCH;
    }

    bool IsMismatch() const noexcept {
        return Operation() == EDLIB_EDOP_MISMATCH;
    }

    bool IsMatchOrMismatch() const noexcept {
        return IsMatch() || IsMismatch();
    }

    // void print() const noexcept {
    //     std::cout << "alignment_length: " << alignment_length << '\n';
    //     std::cout << "align_pos: " << align_pos << '\n';
    //     std::cout << "s1_pos: " << s1_pos << '\n';
    //     std::cout << "s2_pos: " << s2_pos << '\n';
    // }
};

boost::optional<std::string> MakeConsensusString(const std::string & ref, const std::string & a, const std::string & b) {
    AlignResultGuard a_to_b_alignment(EditDistance(a, b, edlib::EDLIB_TASK_PATH));
    if (a_to_b_alignment.status != edlib::EDLIB_STATUS_OK || a_to_b_alignment.editDistance < 0)
        return {};
    // return {a};

    AlignResultGuard a_to_ref_alignment(EditDistance(a, ref, edlib::EDLIB_TASK_PATH));
    AlignResultGuard b_to_ref_alignment(EditDistance(b, ref, edlib::EDLIB_TASK_PATH));
    VERIFY(a_to_ref_alignment.status == edlib::EDLIB_STATUS_OK);
    VERIFY(b_to_ref_alignment.status == edlib::EDLIB_STATUS_OK);

    std::stringstream ss;
    AlignmentIterator a_to_ref(a_to_ref_alignment);
    AlignmentIterator b_to_ref(b_to_ref_alignment);
    AlignmentIterator a_to_b(a_to_b_alignment);

    int previous_unseen_ref_pos = 0;
    bool skipped = false;
    for (; !a_to_b.IsEnd(); ++a_to_b) {
        if (!a_to_b.IsMatchOrMismatch()) {
            skipped = true;
            continue;
        }
        VERIFY(a_to_ref.s1_pos <= a_to_b.s1_pos);
        VERIFY(b_to_ref.s1_pos <= a_to_b.s2_pos);
        while (a_to_ref.s1_pos < a_to_b.s1_pos)
            ++a_to_ref;
        while (b_to_ref.s1_pos < a_to_b.s2_pos)
            ++b_to_ref;
        VERIFY(a_to_ref.s1_pos == a_to_b.s1_pos);
        VERIFY(b_to_ref.s1_pos == a_to_b.s2_pos);

        if (a_to_ref.s2_pos == b_to_ref.s2_pos) {
            auto current_ref_pos = a_to_ref.s2_pos;
            if (skipped && current_ref_pos > previous_unseen_ref_pos) {
                ss << ref.substr(previous_unseen_ref_pos, current_ref_pos - previous_unseen_ref_pos);
            }
            if (a_to_b.IsMatch() || !a_to_ref.IsMatchOrMismatch()) {
                ss << a[a_to_ref.s1_pos];
            } else {
                VERIFY(a_to_b.IsMismatch());
                ss << ref[current_ref_pos];
            }
            previous_unseen_ref_pos = current_ref_pos + 1;
            skipped = false;
        }
    }
    if (previous_unseen_ref_pos != ref.size())
        ss << ref.substr(previous_unseen_ref_pos);

    return { ss.str() };
}

void test() {
    using namespace std;
    static bool call_once = true;
    if (!call_once)
        return;
    call_once = false;

    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAAAGGGGGG"s;
        auto b = "GGGGGGAAAGGGGGG"s;
        auto r = "GGGGGGGGGGGG"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGGGGGGG"s;
        auto b = "GGGGGGGGGGGG"s;
        auto r = "GGGGGGAAAGGGGGG"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAGGGGGG"s;
        auto b = "GGGGGGTGGGGGG"s;
        auto r = "GGGGGGCGGGGGG"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAAAGGGGGG"s;
        auto b = "GGGGGGTTTGGGGGG"s;
        auto r = "GGGGGGGGGGGG"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAAAGGGGGG"s;
        auto b = "GGGGGGAAAGGGGGG"s;
        auto r = "GGGGGGCGGGGGG"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGTGGGGGG"s;
        auto b = "GGGGGGTGGGGGG"s;
        auto r = "GGGGGGAAAGGGGGG"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAAAGGGGGG"s;
        auto b = "GGGGGGTTAGGGGGG"s;
        auto r = "GGGGGGCGGGGGG"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAACGGGGGG"s;
        auto b = "GGGGGGTTCGGGGGG"s;
        auto r = "GGGGGGCGGGGGG"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
    {
        auto a = "AATTCTCTGTGTTGGGGCCCACACCCCAACTTGCATTGCCTGTAGAATTTCTTTTCGAAATTCTCTGTGTTGGGGCCCCTGACTAGAATTGAAAAAAGCTTGTTACAAGCGCATTTTC"s;
        auto b = "AATTCTCTGTGTTGGGGCCCACACCCCAACTTGCATTGTCTGTAGAAATTGGGAATCCAATTTCTCTTTGTTGGGGCCCCTGACTAGAATTGAAAAAAGCTTGTTACAAGCGCATTTTC"s;
        auto r = "AATTCTCTGTGTTGGGGCCCCTGACTAGAATTGAAAAAGCTTGTTACAAGCGCATTTTC"s;
        auto cons = MakeConsensusString(r, a, b);
        VERIFY(cons.is_initialized());
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        std::cout << *cons << '\n';
    }
    std::cout << "=================================\n";
}

} // namespace

boost::optional<std::string> AlignerFiller::operator()(std::vector<Path> const & paths, std::string const & ref) const {
    // test();
    VERIFY(!ref.empty());
    path_extend::ScaffoldSequenceMaker seq_maker(graph);
    std::vector<std::pair<size_t, size_t>> scores; // [(edit distance, path data index)]
    std::vector<std::string> path_data;
    scores.reserve(paths.size());
    path_data.reserve(paths.size());

    #pragma omp parallel for schedule(runtime)
    for (size_t i = 0; i < paths.size(); ++i) {
        auto path = path_extend::BidirectionalPath::create(graph, paths[i]);
        auto query_seq = seq_maker.MakeSequence(*path);
        VERIFY(!query_seq.empty());
        AlignResultGuard edit_distance(EditDistance(query_seq, ref));

        if (edit_distance.status == edlib::EDLIB_STATUS_OK && edit_distance.editDistance >= 0) {
            #pragma omp critical
            {
                scores.emplace_back(edit_distance.editDistance, scores.size());
                path_data.push_back(std::move(query_seq));
            }
        }
    }

    if (scores.empty())
        return {};

    std::sort(scores.begin(), scores.end());

    auto IsGoodScore = [&](auto distance){ return math::le<double, double>(distance, ref.size() * 0.005); };
    if (scores.size() == 1 || IsDominantScore(scores[0].first, scores[1].first)) {
        if (scores.size() > 1)
            TRACE("domination! score[0] = " << scores[0].first << "; score[1] = " << scores[1].first << "; amount_of_paths = " << scores.size());
        if (IsGoodScore(scores.front().first))
            return {std::move(path_data[scores.front().second])};
        return {};
    }

    TRACE("domination failed! score[0] = " << scores[0].first << "; score[1] = " << scores[1].first << "; amount_of_paths = " << scores.size());
    if (scores.size() == 2) {
        static std::ofstream out1("best_of_two");
        static std::ofstream out2("worst_of_two");
        static std::ofstream out_ref("ref_seq");
        static std::ofstream out_cons("consensus");
        static size_t id = 0;
        TRACE("Current id = " << id);

        out1 << ">case" << id << '\n';
        out1 << path_data[scores[0].second] << '\n';
        out2 << ">case" << id << '\n';
        out2 << path_data[scores[1].second] << '\n';
        out_ref << ">case" << id << '\n';
        out_ref << ref << '\n';
        out_cons << ">case" << id << '\n';
        ++id;

        TRACE("Trying to make consensus");
        // return MakeConsensusString(ref, scores[0], scores[1]);
        auto cons = MakeConsensusString(ref, path_data[scores[0].second], path_data[scores[1].second]);
        if (cons.is_initialized()) {
            AlignResultGuard res(EditDistance(*cons, ref));
            TRACE("Succsess! New score: " << res.editDistance);
            out_cons << *cons << '\n';
            if (!IsGoodScore(res.editDistance)) {
                TRACE("But the distance is too high!");
                cons.reset();
                VERIFY(!cons.is_initialized());
            }
        } else {
            TRACE("Failed!");
            out_cons << "NOTHING" << '\n';
        }
        // return cons;
        return {};
    }
    return {};
}

bool AlignerFiller::IsDominantScore(size_t dominator, size_t other) const noexcept {
    return dominator + 10 < other &&
            math::ls((double)dominator, (double)other * score_domination_coeff);
}

} // namespace sequence_corrector
