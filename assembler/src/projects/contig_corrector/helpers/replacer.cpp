#include "replacer.hpp"
#include "string_utils.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"
#include "helpers/replace_info_writer.hpp"

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
    VERIFY(contig_start_pos <= contig_end_pos);
}

void ReplaceInfo::DropBehind(size_t size) noexcept {
    drop_from_tail += size;
    contig_end_pos -= size;
    VERIFY(contig_start_pos <= contig_end_pos);
}

namespace {

void DumpReplaceInfo(string const & seq_name, string const & seq, list<ReplaceInfo> const & mapping_info) {
    if (!ReplaceInfoWriter::stream)
        return;

    auto & stream = *ReplaceInfoWriter::stream;
    auto Skip = [](size_t & i, auto const & s, int (*pred)(int)) {
        while (i + s.drop_from_tail < s.seq.size() && pred(s.seq[i]))
            ++i;
    };
    size_t current_pos = 0;
    size_t result_pos = 0;
    for (auto const & path : mapping_info) {
        {
            auto len = path.contig_start_pos - current_pos;
            if (len > 0) {
                stream.Write(seq_name, RangeType::origin, result_pos, len);
                result_pos += len;
            }
        }

        for (size_t i = path.drop_from_head; i + path.drop_from_tail < path.seq.size();) {
            auto start_pos = i;
            if (islower(path.seq[i])) {
                Skip(i, path, islower);
                stream.Write(seq_name, RangeType::edge, result_pos, i - start_pos);
            } else {
                Skip(i, path, isupper);
                stream.Write(seq_name, RangeType::path, result_pos, i - start_pos);
            }
            result_pos += i - start_pos;
        }

        current_pos = path.contig_end_pos;
    }
    if (current_pos < seq.size()) {
        auto len = seq.size() - current_pos;
        stream.Write(seq_name, RangeType::origin, result_pos, len);
    }
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
    // we use a strategy when overlapped parts will be dropped from both paths
    size_t tail_drop = 0;
    while (true) {
        auto next_path = std::next(current_path);
        if (next_path == mapping_info.end())
            break;
        VERIFY_MSG(current_path->contig_start_pos <= next_path->contig_start_pos, "This should not have happened! :(");
        if (next_path->contig_start_pos < current_path->contig_end_pos) {
            auto overlapping = current_path->contig_end_pos - next_path->contig_start_pos;
            tail_drop += overlapping; // this should be max function!
            next_path->DropAhead(overlapping);
        }
        if (next_path->ShouldBeDropped()) {
            VERIFY_MSG(false, "Wrong behavior");
            mapping_info.erase(next_path);
            continue;
        }

        current_path->DropBehind(tail_drop);
        tail_drop = 0;
        if (current_path->ShouldBeDropped()) {
            VERIFY_MSG(false, "This should not have happened! :(");
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
