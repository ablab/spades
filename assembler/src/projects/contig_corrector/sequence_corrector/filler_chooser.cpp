#include "filler_chooser.hpp"
#include "common/assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include "utils/logger/logger.hpp"

#include "edlib/edlib.h"

#include <utility>
#include <algorithm>
#include <sstream>
#include <unordered_set>

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

    bool IsInsertedToFirst() const noexcept {
        return Operation() == EDLIB_EDOP_INSERT;
    }

    bool IsInsertedToSecond() const noexcept {
        return Operation() == EDLIB_EDOP_DELETE;
    }

    // void print() const noexcept {
    //     std::cout << "alignment_length: " << alignment_length << '\n';
    //     std::cout << "align_pos: " << align_pos << '\n';
    //     std::cout << "s1_pos: " << s1_pos << '\n';
    //     std::cout << "s2_pos: " << s2_pos << '\n';
    // }
};

void SkipNotMatched(AlignmentIterator & it) {
    while (!it.IsEnd() && !it.IsMatch())
        ++it;
}

boost::optional<std::string> MakeCommonStr(const std::string & a, const std::string & b) {
    AlignResultGuard a_to_b_alignment(EditDistance(a, b, edlib::EDLIB_TASK_PATH));
    if (a_to_b_alignment.status != edlib::EDLIB_STATUS_OK || a_to_b_alignment.editDistance < 0)
        return {};

    std::stringstream ss;
    AlignmentIterator a_to_b(a_to_b_alignment);
    while (!a_to_b.IsEnd()) {
        if (a_to_b.IsMatch()) {
            ss << a[a_to_b.s1_pos];
            ++a_to_b;
        } else {
            ss << '.';
            SkipNotMatched(a_to_b);
        }
    }

    return ss.str();
}

boost::optional<std::string> MakeConsensusString(const std::string & ref, const std::string & a, const std::string & b) {
    auto common_str = MakeCommonStr(a, b);
    if (!common_str) {
        WARN("a/b distance too large");
        return {};
    }

    const std::string & c = *common_str;

    AlignResultGuard c_to_ref_alignment(EditDistance(c, ref, edlib::EDLIB_TASK_PATH));
    if (c_to_ref_alignment.status != edlib::EDLIB_STATUS_OK || c_to_ref_alignment.editDistance < 0) {
        WARN("common_str/ref distance too large");
        return {};
    }

    std::stringstream ss;
    AlignmentIterator c_to_ref(c_to_ref_alignment);
    while (!c_to_ref.IsEnd()) {
        if (c[c_to_ref.s1_pos] != '.') {
            if (!c_to_ref.IsInsertedToSecond())
                ss << c[c_to_ref.s1_pos];
            ++c_to_ref;
            continue;
        }

        for (;!c_to_ref.IsEnd() && c[c_to_ref.s1_pos] == '.'; ++c_to_ref) {
            if (!c_to_ref.IsInsertedToFirst())
                ss << ref[c_to_ref.s2_pos];
        }
        for (;!c_to_ref.IsEnd() && c_to_ref.IsInsertedToSecond(); ++c_to_ref)
            ss << ref[c_to_ref.s2_pos];
    }

    return ss.str();
}

boost::optional<std::string> GetConsStrForTwoPaths(const std::string & ref, const std::string & a, const std::string & b, const std::function<bool(int)> & IsGoodScore) {
    static std::ofstream out1("best_of_two");
    static std::ofstream out2("worst_of_two");
    static std::ofstream out_ref("ref_seq");
    static std::ofstream out_cons("consensus");
    static size_t id = 0;
    WARN("Current id = " << id);

    out1 << ">case" << id << '\n';
    out1 << a << '\n';
    out2 << ">case" << id << '\n';
    out2 << b << '\n';
    out_ref << ">case" << id << '\n';
    out_ref << ref << '\n';
    out_cons << ">case" << id << '\n';
    ++id;

    WARN("Trying to make consensus");
    auto cons = MakeConsensusString(ref, a, b);
    if (cons.is_initialized()) {
        AlignResultGuard res(EditDistance(*cons, ref));
        WARN("Succsess! New score: " << res.editDistance);
        out_cons << *cons << '\n';
        if (!IsGoodScore(res.editDistance)) {
            WARN("But the distance is too high!");
            cons.reset();
        }
    } else {
        WARN("Failed!");
        out_cons << "NOTHING" << '\n';
    }
    return cons;
}

} // namespace

struct Node {
    sequence_corrector::Path::value_type edge;
    size_t id;

    bool operator == (const Node& other) const noexcept {
        return edge == other.edge && id == other.id;
    }
};

std::ostream& operator << (std::ostream& out, const Node& node) {
    out << node.edge << '/' << node.id;
    return out;
}

namespace std {

template<>
struct hash<Node> {
    size_t operator()(const Node& node) const noexcept {
        size_t seed = hash<decltype(node.edge)>()(node.edge);
        seed ^= hash<decltype(node.id)>()(node.id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

} // namespace std

namespace sequence_corrector {

namespace {

struct TwoEdgeSet {
    boost::optional<Node> first_edge;
    boost::optional<Node> second_edge;

    bool TryInsert(size_t id, Path::value_type edge) noexcept {
        auto TrySet = [&](boost::optional<Node> & this_edge) noexcept {
            if (!this_edge) {
                this_edge = Node{edge, id};
                return true;
            }

            return this_edge->edge == edge && this_edge->id == id;
        };

        return TrySet(first_edge) || TrySet(second_edge);
    }
};

std::ostream& operator << (std::ostream& out, const TwoEdgeSet& set) {
    out << "{ ";

    if (set.first_edge)
        out << *set.first_edge;
    else
        out << "??";

    out << " | ";

    if (set.second_edge)
        out << *set.second_edge;
    else
        out << "??";
    out << " }";
    return out;
}

bool Find(std::vector<Path> const & paths, Path const & path) {
    for (auto const & sample : paths) {
        if (sample == path)
            return true;
    }
    return false;
}

using PathGraph = std::unordered_map<Node, TwoEdgeSet>;

boost::optional<std::pair<Path, Path>> TryReversedCompress(std::vector<Path> const & paths);

boost::optional<std::pair<Path, Path>> Compress(std::vector<Path> const & paths, bool reversed = false) {
    PathGraph graph;

    auto DumpPaths = [&](auto const & msg) {
        std::cout << msg << " (paths size: " << paths.size() << ") :\n";
        if (paths.size() > 16)
            return;
        std::cout << "-------------------\n";
        for (auto const & path : paths) {
            std::unordered_map<Path::value_type, size_t> free_id;
            for (auto const & edge : path)
                std::cout << edge << '/' << free_id[edge]++ << ' ';
            std::cout << '\n';
        }
        std::cout << "-------------------\n";
        std::cout << "Graph:\n";
        for (auto const & el : graph)
            std::cout << el.first << " -> " << el.second << '\n';
        std::cout << "-------------------\n";
    };

    std::unordered_set<Node> unused_nodes;
    TwoEdgeSet entry_points;
    for (auto const & path : paths) {
        if (!entry_points.TryInsert(0, path.front()))
            return DumpPaths("Too many roots"), boost::none;
        unused_nodes.insert({path.front(), 0});
    }

    for (auto const & path : paths) {
        std::unordered_map<Path::value_type, size_t> free_id;
        Node prev_edge = {path.front(), free_id[path.front()]++};
        for (size_t i = 1; i < path.size(); ++i) {
            auto edge = path[i];
            auto id = free_id[edge]++;
            if (!graph[prev_edge].TryInsert(id, edge))
                return DumpPaths("Too many ways, " + 
                                 std::to_string(prev_edge.edge.id_) + '/' + std::to_string(prev_edge.id) +
                                 " -?-> " + std::to_string(edge.id_) + '/' + std::to_string(id)), (reversed ? boost::none : TryReversedCompress(paths));
            prev_edge = {edge, id};
            unused_nodes.insert(prev_edge);
        }
    }

    std::pair<Path, Path> ans;
    std::unordered_map<Node, unsigned char> used;
    Node* cur_edge = &*entry_points.first_edge;
    while (true) {
        unused_nodes.erase(*cur_edge);
        auto& cur_used_cnt = ++used[*cur_edge];
        if (cur_used_cnt > 1)
            return DumpPaths("Loop is found"), (reversed ? boost::none : TryReversedCompress(paths)); // loop is found
        ans.first.push_back(cur_edge->edge);
        auto it = graph.find(*cur_edge);
        if (it == graph.end())
            break;
        cur_edge = &*it->second.first_edge;
    }

    cur_edge = (entry_points.second_edge ? &*entry_points.second_edge : &*entry_points.first_edge);
    while (true) {
        unused_nodes.erase(*cur_edge);
        auto& cur_used_cnt = ++used[*cur_edge];
        if (cur_used_cnt > 2)
            return DumpPaths("Loop is found"), (reversed ? boost::none : TryReversedCompress(paths)); // loop is found
        ans.second.push_back(cur_edge->edge);
        auto it = graph.find(*cur_edge);
        if (it == graph.end())
            break;
        if (!it->second.second_edge) {
            cur_edge = &*it->second.first_edge;
            continue;
        }
        if (used[*it->second.first_edge] == 0)
            return DumpPaths("Cannot compress to 2 paths"), (reversed ? boost::none : TryReversedCompress(paths));
        cur_edge = &*it->second.second_edge;
    }

    if (!Find(paths, ans.first) || !Find(paths, ans.second))
        return DumpPaths("Too strange paths"), (reversed ? boost::none : TryReversedCompress(paths));

    if (reversed) {
        if (!unused_nodes.empty())
            std::cout << "There are unused nodes!\n";
        DumpPaths("Reversed graph!");
    }

    if (!unused_nodes.empty())
        return boost::none;

    return {std::move(ans)};
}

/// TODO: should we really do that?
boost::optional<std::pair<Path, Path>> TryReversedCompress(std::vector<Path> const & paths) {
    std::vector<Path> reversed_paths(paths.size());
    for (size_t i = 0; i < paths.size(); ++i)
        reversed_paths[i].insert(reversed_paths[i].end(), paths[i].crbegin(), paths[i].crend());
    std::cout << "Trying to compress via reverse\n";
    auto compressed_paths = Compress(reversed_paths, true);
    if (compressed_paths) {
        std::cout << "Success reversed compressing!\n";
        std::reverse(compressed_paths->first.begin(), compressed_paths->first.end());
        std::reverse(compressed_paths->second.begin(), compressed_paths->second.end());
    } else {
        std::cout << "Failed reversed compressing!\n";
    }
    return compressed_paths;
}


struct Score {
    int distance;
    size_t path_data_index;
    size_t origin_path_index;

    bool operator < (const Score& other) const noexcept {
        if (distance != other.distance)
            return distance < other.distance;
        if (path_data_index != other.path_data_index)
            return path_data_index < other.path_data_index;
        return origin_path_index < other.origin_path_index;
    }
};

std::vector<Path> Filter(const std::vector<Path>& paths, const std::vector<Score>& scores) {
    std::vector<Path> filtered_paths;
    filtered_paths.reserve(scores.size());
    for (const auto& path : scores)
        filtered_paths.push_back(paths[path.origin_path_index]);
    return filtered_paths;
}

} // namespace

void CompressTest() {
    auto Check = [](const std::vector<Path>& paths) {
        std::cout << "paths:\n";
        for (const auto& path : paths) {
            for (auto const & edge : path)
                std::cout << edge << ' ';
            std::cout << '\n';
        }
        std::cout.flush();
        auto compressed = Compress(paths);
        if (!compressed) {
            std::cout << "Cannot compress!\n";
            return;
        }
        std::cout << "compressed:\n";
        for (auto const & edge : compressed->first)
            std::cout << edge << ' ';
        std::cout << '\n';
        for (auto const & edge : compressed->second)
            std::cout << edge << ' ';
        std::cout << "\n--------------------------\n";
    };

    Check({{1,2,3}, {1,2,3}});
    Check({{1,2,3}, {1,4,3}});
    Check({{1,2,3}, {1,3,3}});
    Check({{1}, {2}});
    Check({{1,3}, {2,4}});
    Check({{1,5,3}, {2,5,4}});
    Check({{1,10,2,11,3}, {1,12,2,11,3}, {1,10,2,13,3}, {1,12,2,13,3}});
}

void ConsensusTest() {
    using namespace std;

    auto Check = [](const auto & r, const auto & a, const auto & b) {
        auto cons = MakeConsensusString(r, a, b);
        std::cout << a << '\n';
        std::cout << b << '\n';
        std::cout << r << '\n';
        if (cons)
            std::cout << *cons << '\n';
        else
            std::cout << "no consensus provided\n";
    };

    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAAAGGGGGG"s;
        auto b = "GGGGGGAAAGGGGGG"s;
        auto r = "GGGGGGGGGGGG"s;
        Check(r, a, b);

    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGGGGGGG"s;
        auto b = "GGGGGGGGGGGG"s;
        auto r = "GGGGGGAAAGGGGGG"s;
        Check(r, a, b);
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAGGGGGG"s;
        auto b = "GGGGGGTGGGGGG"s;
        auto r = "GGGGGGCGGGGGG"s;
        Check(r, a, b);
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAAAGGGGGG"s;
        auto b = "GGGGGGTTTGGGGGG"s;
        auto r = "GGGGGGGGGGGG"s;
        Check(r, a, b);
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAAAGGGGGG"s;
        auto b = "GGGGGGAAAGGGGGG"s;
        auto r = "GGGGGGCGGGGGG"s;
        Check(r, a, b);
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGTGGGGGG"s;
        auto b = "GGGGGGTGGGGGG"s;
        auto r = "GGGGGGAAAGGGGGG"s;
        Check(r, a, b);
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAAAGGGGGG"s;
        auto b = "GGGGGGTTAGGGGGG"s;
        auto r = "GGGGGGCGGGGGG"s;
        Check(r, a, b);
    }
    std::cout << "=================================\n";
    {
        auto a = "GGGGGGAACGGGGGG"s;
        auto b = "GGGGGGTTCGGGGGG"s;
        auto r = "GGGGGGCGGGGGG"s;
        Check(r, a, b);
    }
    std::cout << "=================================\n";
    {
        auto a = "AATTCTCTGTGTTGGGGCCCACACCCCAACTTGCATTGCCTGTAGAATTTCTTTTCGAAATTCTCTGTGTTGGGGCCCCTGACTAGAATTGAAAAAAGCTTGTTACAAGCGCATTTTC"s;
        auto b = "AATTCTCTGTGTTGGGGCCCACACCCCAACTTGCATTGTCTGTAGAAATTGGGAATCCAATTTCTCTTTGTTGGGGCCCCTGACTAGAATTGAAAAAAGCTTGTTACAAGCGCATTTTC"s;
        auto r = "AATTCTCTGTGTTGGGGCCCCTGACTAGAATTGAAAAAGCTTGTTACAAGCGCATTTTC"s;
        Check(r, a, b);
    }
    std::cout << "=================================\n";
}

boost::optional<std::string> AlignerFiller::operator()(std::vector<Path> const & paths, std::string const & ref) const {
    VERIFY(!ref.empty());
    path_extend::ScaffoldSequenceMaker seq_maker(graph);
    std::vector<Score> scores; // [(edit distance, path data index)]
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
                scores.push_back({edit_distance.editDistance, scores.size(), i});
                path_data.push_back(std::move(query_seq));
            }
        }
    }

    std::cout << "ref.size: " << ref.size() << ", scores.size: " << scores.size() << ", paths.size: " << paths.size() << '\n';

    if (scores.empty())
        return boost::none;

    std::sort(scores.begin(), scores.end());

    if (3 <= scores.size() && scores.size() <= 9) {
        for (auto const & path : scores) {
            std::cout << "[" << path.distance << "] ";
            for (auto const & edge : paths[path.path_data_index])
                std::cout << edge << "{" << graph.length(edge) << "}" << ' ';
            std::cout << '\n';
        }
    }

    auto IsGoodScore = [&](auto distance){ return math::le<double, double>(distance, ref.size() * 0.02); };
    if (scores.size() == 1 || IsDominantScore(scores[0].distance, scores[1].distance)) {
        if (scores.size() == 1)
            std::cout << "the score is " << scores[0].distance << "\n";
        if (scores.size() > 1)
            std::cout << "domination! score[0] = " << scores[0].distance << "; score[1] = " << scores[1].distance << "; amount_of_paths = " << scores.size() << '\n';
        if (IsGoodScore(scores.front().distance))
            return {std::move(path_data[scores.front().path_data_index])};
        std::cout << "But the score is bad\n";
        return boost::none;
    }

    std::cout << "domination failed! score[0] = " << scores[0].distance << "; score[1] = " << scores[1].distance << "; amount_of_paths = " << scores.size() << '\n';

    if (scores.size() == 2)
        return GetConsStrForTwoPaths(ref, path_data[scores[0].path_data_index], path_data[scores[1].path_data_index], IsGoodScore);

    std::cout << "trying to compress " << scores.size() << " paths" << std::endl;
    boost::optional<std::pair<Path, Path>> compressed_paths;
    if (scores.size() == paths.size())
        compressed_paths = Compress(paths);
    else
        compressed_paths = Compress(Filter(paths, scores));

    if (!compressed_paths) {
        std::cout << "fail!" << std::endl;
        return boost::none;
    }

    std::cout << "compressed paths:\n";
    for (auto const & edge : compressed_paths->first)
        std::cout << edge << "{" << graph.length(edge) << "}" << ' ';
    std::cout << '\n';
    for (auto const & edge : compressed_paths->second)
        std::cout << edge << "{" << graph.length(edge) << "}" << ' ';
    std::cout << '\n';

    auto path1 = path_extend::BidirectionalPath::create(graph, compressed_paths->first);
    auto query_seq1 = seq_maker.MakeSequence(*path1);
    auto path2 = path_extend::BidirectionalPath::create(graph, compressed_paths->second);
    auto query_seq2 = seq_maker.MakeSequence(*path2);

    auto big_cons = GetConsStrForTwoPaths(ref, query_seq1, query_seq2, IsGoodScore);
    if (big_cons)
        std::cout << "success!" << std::endl;

    return big_cons;
}

bool AlignerFiller::IsDominantScore(size_t dominator, size_t other) const noexcept {
    return dominator + 10 < other &&
            math::ls((double)dominator, (double)other * score_domination_coeff);
}

} // namespace sequence_corrector
