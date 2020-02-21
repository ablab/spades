#include "kmer_map.hpp"

namespace debruijn_graph {

std::size_t KMerMap::str_hash::operator()(const char* key, std::size_t key_size) const {
    return CityHash64(key, key_size);
}

KMerMap::iterator::iterator(unsigned k, HTMap::const_iterator iter)
    : k_(k)
    , iter_(iter)
{}

void KMerMap::iterator::increment() {
    ++iter_;
}

bool KMerMap::iterator::equal(const iterator &other) const {
    return iter_ == other.iter_;
}

auto KMerMap::iterator::dereference() const -> const std::pair<Kmer, Seq> {
    iter_.key(key_out_);
    Kmer k(k_, (const RawSeqData*)key_out_.data());
    Seq s(k_, (const RawSeqData*)iter_.value());
    return std::make_pair(k, s);
}

KMerMap::KMerMap(unsigned k)
    : k_(k)
{
    rawcnt_ = (unsigned)Seq::GetDataSize(k_);
}

KMerMap::~KMerMap() {
    clear();
}

void KMerMap::erase(const Kmer &key) {
    auto res = mapping_.find_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData));
    if (res == mapping_.end())
        return;

    delete[] res.value();
    mapping_.erase(res);
}

void KMerMap::set(const Kmer &key, const Seq &value) {
    RawSeqData *rawvalue = nullptr;
    auto res = mapping_.find_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData));
    if (res == mapping_.end()) {
        rawvalue = new RawSeqData[rawcnt_];
        mapping_.insert_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData), rawvalue);
    } else {
        rawvalue = res.value();
    }
    memcpy(rawvalue, value.data(), rawcnt_ * sizeof(RawSeqData));
}

bool KMerMap::count(const Kmer &key) const {
    return mapping_.count_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData));
}

auto KMerMap::find(const Kmer &key) const -> const RawSeqData* {
    auto res = mapping_.find_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData));
    if (res == mapping_.end())
        return nullptr;

    return res.value();
}

auto KMerMap::find(const RawSeqData *key) const -> const RawSeqData* {
    auto res = mapping_.find_ks((const char*)key, rawcnt_ * sizeof(RawSeqData));
    if (res == mapping_.end())
        return nullptr;

    return res.value();
}

void KMerMap::clear() {
    // Delete all the values
    for (auto it = mapping_.begin(); it != mapping_.end(); ++it) {
        VERIFY(it.value() != nullptr);
        delete[] it.value();
        it.value() = nullptr;
    }

    // Delete the mapping and all the keys
    mapping_.clear();
}

size_t KMerMap::size() const {
    return mapping_.size();
}

auto KMerMap::begin() const -> iterator {
    return iterator(k_, mapping_.begin());
}

auto KMerMap::end() const -> iterator {
    return iterator(k_, mapping_.end());
}

} // namespace debruijn_graph
