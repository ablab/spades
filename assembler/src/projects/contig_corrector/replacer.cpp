#include "replacer.hpp"
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <sstream>

using namespace std;

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

string MakeSeq(string const & seq, list<ReplaceInfo> const & mapping_info) {
    size_t current_pos = 0;
    stringstream ss;
    size_t i = 0;
    for (auto const & path : mapping_info) {
        if (path.drop_from_head > 0)
            WARN("The new path prefix with len=" << path.drop_from_head << " of " << path.Size() << " would be dropped");
        if (path.drop_from_tail > 0)
            WARN("The new path suffix with len=" << path.drop_from_tail << " of " << path.Size() << " would be dropped");
        ss << seq.substr(current_pos, path.contig_start_pos - current_pos);
        ss << path.seq.substr(path.drop_from_head, path.ResultSize());
        current_pos = path.contig_end_pos;
    }
    ss << seq.substr(current_pos);
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
