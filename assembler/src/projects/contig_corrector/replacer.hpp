#pragma once
#include <string>
#include <list>

struct ReplaceInfo {
    std::string seq;
    size_t contig_start_pos; // inclusive
    size_t contig_end_pos;   // exclusive
    size_t drop_from_head;
    size_t drop_from_tail;

    ReplaceInfo(std::string seq, size_t contig_start_pos = 0, size_t contig_end_pos = 0, size_t drop_from_head = 0, size_t drop_from_tail = 0)
        : seq(std::move(seq))
        , contig_start_pos(contig_start_pos)
        , contig_end_pos(contig_end_pos)
        , drop_from_head(drop_from_head)
        , drop_from_tail(drop_from_tail)
    {}

    size_t Size() const noexcept { return seq.size(); }

    size_t StartPos() const noexcept { return contig_start_pos + drop_from_head; }
    size_t EndPos() const noexcept { return contig_end_pos - drop_from_tail; }

    bool ShouldBeDropped() const noexcept {
        return Size() <= drop_from_head + drop_from_tail;
    }

    size_t ResultSize() const noexcept;

    void DropAhead(size_t size) noexcept;
    void DropBehind(size_t size) noexcept;
};

std::string Replace(std::string const & seq, std::list<ReplaceInfo> && replace_info);
