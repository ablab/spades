#pragma once
#include <string>
#include <list>
#include "utils/verify.hpp"

namespace helpers {

struct ReplaceInfo {
    ///@invariant \ref seq [ \ref drop_from_head ] points at \ref contig_start_pos in original sequence
    ///@invariant \ref seq [ seq.size() - \ref drop_from_tail - 1 ] points at \ref contig_end_pos - 1 in original sequence

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
    {
        VERIFY(contig_start_pos <= contig_end_pos);
    }

    size_t Size() const noexcept {
        VERIFY(contig_start_pos <= contig_end_pos);
        return seq.size();
    }

    bool ShouldBeDropped() const noexcept {
        VERIFY(contig_start_pos <= contig_end_pos);
        return Size() <= drop_from_head + drop_from_tail;
    }

    size_t ResultSize() const noexcept;

    void DropAhead(size_t size) noexcept;
    void DropBehind(size_t size) noexcept;
};

std::string Replace(std::string const & seq, std::list<ReplaceInfo> && replace_info);
std::string ReplaceAndDump(std::string const & seq, std::list<ReplaceInfo> && replace_info, std::string const & seq_name);

} // namespace helpers
