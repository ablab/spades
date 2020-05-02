//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence_mapper_notifier.hpp"

#include "sequence_mapper.hpp"
#include "io/reads/paired_read.hpp"
#include "io/reads/read_stream_vector.hpp"

namespace debruijn_graph {

SequenceMapperNotifier::SequenceMapperNotifier(const GraphPack& gp, size_t lib_count)
    : gp_(gp)
    , listeners_(lib_count) 
{}

void SequenceMapperNotifier::Subscribe(size_t lib_index, SequenceMapperListener* listener) {
    VERIFY(lib_index < listeners_.size());
    listeners_[lib_index].push_back(listener);
}

void SequenceMapperNotifier::NotifyStartProcessLibrary(size_t ilib, size_t thread_count) const {
    for (const auto& listener : listeners_[ilib])
        listener->StartProcessLibrary(thread_count);
}

void SequenceMapperNotifier::NotifyStopProcessLibrary(size_t ilib) const {
    for (const auto& listener : listeners_[ilib])
        listener->StopProcessLibrary();
}

void SequenceMapperNotifier::NotifyMergeBuffer(size_t ilib, size_t ithread) const {
    std::string thread_str = std::to_string(ithread);
    TIME_TRACE_SCOPE("SequenceMapperNotifier::MergeBuffer", thread_str);
    for (const auto& listener : listeners_[ilib])
        listener->MergeBuffer(ithread);
}

template<>
void SequenceMapperNotifier::NotifyProcessRead(const io::PairedReadSeq& r,
                                               const SequenceMapperT& mapper,
                                               size_t ilib,
                                               size_t ithread) const
{
    const Sequence& read1 = r.first().sequence();
    const Sequence& read2 = r.second().sequence();
    MappingPath<EdgeId> path1 = mapper.MapSequence(read1);
    MappingPath<EdgeId> path2 = mapper.MapSequence(read2);
    for (const auto& listener : listeners_[ilib]) {
        listener->ProcessPairedRead(ithread, r, path1, path2);
        listener->ProcessSingleRead(ithread, r.first(), path1);
        listener->ProcessSingleRead(ithread, r.second(), path2);
    }
}

template<>
void SequenceMapperNotifier::NotifyProcessRead(const io::PairedRead& r,
                                               const SequenceMapperT& mapper,
                                               size_t ilib,
                                               size_t ithread) const
{
    MappingPath<EdgeId> path1 = mapper.MapRead(r.first());
    MappingPath<EdgeId> path2 = mapper.MapRead(r.second());
    for (const auto& listener : listeners_[ilib]) {
        listener->ProcessPairedRead(ithread, r, path1, path2);
        listener->ProcessSingleRead(ithread, r.first(), path1);
        listener->ProcessSingleRead(ithread, r.second(), path2);
    }
}

template<>
void SequenceMapperNotifier::NotifyProcessRead(const io::SingleReadSeq& r,
                                               const SequenceMapperT& mapper,
                                               size_t ilib,
                                               size_t ithread) const
{
    const Sequence& read = r.sequence();
    MappingPath<EdgeId> path = mapper.MapSequence(read);
    for (const auto& listener : listeners_[ilib])
        listener->ProcessSingleRead(ithread, r, path);
}

template<>
void SequenceMapperNotifier::NotifyProcessRead(const io::SingleRead& r,
                                               const SequenceMapperT& mapper,
                                               size_t ilib,
                                               size_t ithread) const
{
    MappingPath<EdgeId> path = mapper.MapRead(r);
    for (const auto& listener : listeners_[ilib])
        listener->ProcessSingleRead(ithread, r, path);
}

} // namespace debruijn_graph
