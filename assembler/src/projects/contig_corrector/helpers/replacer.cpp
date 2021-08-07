#include "replacer.hpp"
#include "string_utils.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <sstream>
#include <iomanip>
#include <fstream>

using namespace std;

namespace helpers {

size_t ReplaceInfo::ResultSize() const noexcept {
    VERIFY_MSG(!ShouldBeDropped(), "LoL? LengthInNucls() = " << Size() << ", drop_from_head = " << drop_from_head << ", drop_from_tail = " << drop_from_tail);
    return Size() - drop_from_head - drop_from_tail;
}

void ReplaceInfo::DropAhead(size_t size) noexcept {
    drop_from_head += size;
    contig_start_pos += size;
}

void ReplaceInfo::DropBehind(size_t size) noexcept {
    drop_from_tail += size;
    contig_end_pos -= size;
}

namespace {

void DumpReplaceInfo(string const & seq_name, string const & seq, list<ReplaceInfo> const & mapping_info) {
    static ofstream out("replace_info_dump.dump");
    static bool once_flag = true;
    auto Print = [&out](auto const & name, auto const & type, auto const & from, auto const & len) {
        out << setw(15) << name << ' ' << setw(10) << type << ' ' << setw(10) << from << ' ' << setw(10) << len << '\n';
    };
    auto Skip = [](size_t & i, auto const & s, int (*pred)(int)) {
        while (i < s.size() && pred(s[i]))
            ++i;
    };
    if (once_flag) {
        Print("name", "seq_type", "from", "len");
        once_flag = false;
    }
    size_t current_pos = 0;
    size_t result_pos = 0;
    size_t total_origin_len = 0;
    size_t total_replaced_len = 0;
    for (auto const & path : mapping_info) {
        {
            auto len = path.contig_start_pos - current_pos;
            if (len > 0) {
                Print(seq_name, "origin", result_pos, len);
                result_pos += len;
                total_origin_len += len;
            }
        }

        for (size_t i = path.drop_from_head; i + path.drop_from_tail < path.seq.size();) {
            auto start_pos = i;
            if (islower(path.seq[i])) {
                Skip(i, path.seq, islower);
                Print(seq_name, "edge", start_pos + result_pos, i - start_pos);
            } else {
                Skip(i, path.seq, isupper);
                Print(seq_name, "path", start_pos + result_pos, i - start_pos);
            }
        }
        auto len = path.ResultSize();
        result_pos += len;
        total_replaced_len += len;

        current_pos = path.contig_end_pos;
    }
    if (current_pos < seq.size()) {
        auto len = seq.size() - current_pos;
        Print(seq_name, "origin", result_pos, len);
       total_origin_len += len;
    }
    INFO("Corrected " << total_replaced_len << " nucls (" << 100.0*total_replaced_len / (total_replaced_len + total_origin_len) << "%) of " << seq_name);
}

string MakeSeq(string const & seq, list<ReplaceInfo> const & mapping_info) {
    size_t current_pos = 0;
    stringstream ss;
    for (auto const & path : mapping_info) {
        if (path.drop_from_head > 0)
            WARN("The new path prefix with len=" << path.drop_from_head << " of " << path.Size() << " would be dropped");
        if (path.drop_from_tail > 0)
            WARN("The new path suffix with len=" << path.drop_from_tail << " of " << path.Size() << " would be dropped");
        ss << ToLowerStr(seq.substr(current_pos, path.contig_start_pos - current_pos));
        ss << ToUpperStr(path.seq.substr(path.drop_from_head, path.ResultSize()));
        current_pos = path.contig_end_pos;
    }
    ss << ToLowerStr(seq.substr(current_pos));
    return ss.str();
}

void DropOverlappedPrefixes(list<ReplaceInfo> & mapping_info) {
    size_t start_size = mapping_info.size();
    if (mapping_info.empty())
        return;

    auto current_path = mapping_info.begin();
    size_t tail_drop = 0;
    while (true) {
        auto next_path = std::next(current_path);
        if (next_path == mapping_info.end())
            break;
        if (next_path->contig_start_pos < current_path->contig_end_pos) {
            auto overlapping = current_path->contig_end_pos - next_path->contig_start_pos;
            tail_drop += overlapping;
            next_path->DropAhead(overlapping);
        }
        if (next_path->ShouldBeDropped()) {
            mapping_info.erase(next_path);
            continue;
        }

        current_path->DropBehind(tail_drop);
        tail_drop = 0;
        if (current_path->ShouldBeDropped()) {
            current_path = mapping_info.erase(current_path);
            continue;
        }
        ++current_path;
    }
    if (mapping_info.size() != start_size) {
        size_t dropped_pats_cnt = start_size - mapping_info.size();
        WARN("Oh no! Full " << dropped_pats_cnt << " path" << (dropped_pats_cnt != 1 ? "s" : "") << " would be dropped");
    }
}

} // namespace

std::string Replace(std::string const & seq, std::list<ReplaceInfo> && replace_info) {
    DropOverlappedPrefixes(replace_info);
    return MakeSeq(seq, replace_info);
}

std::string ReplaceAndDump(std::string const & seq, std::list<ReplaceInfo> && replace_info, std::string const & seq_name) {
    DropOverlappedPrefixes(replace_info);
    DumpReplaceInfo(seq_name, seq, replace_info);
    return MakeSeq(seq, replace_info);
}

} // namespace helpers
