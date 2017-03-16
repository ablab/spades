//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence/nucl.hpp"
#include "io/reads/modifying_reader_wrapper.hpp"
#include "adt/bag.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
//#include "indices/edge_index_builders.hpp"

namespace cap {

struct Count {
    size_t count;
    Count() : count(0) {}
    Count conjugate(size_t /*k*/) {
        return *this;
    }
};

template <class Builder>
class RepeatSearchingIndexBuilder : public Builder {
    typedef Builder base;
 public:
    typedef typename Builder::IndexT IndexT;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::KMerIdx KmerIdx;
    typedef typename IndexT::KeyWithHash KWH;

 private:
    template<class ReadStream>
    size_t FillCoverageFromStream(ReadStream &stream,
                                  IndexT &index) const {
        unsigned k = index.k();
        while (!stream.eof()) {
            typename ReadStream::ReadT r;
            stream >> r;

            const Sequence &seq = r.sequence();
            if (seq.size() < k)
                continue;

            KWH kwh = index.ConstructKWH(seq.start<Kmer>(k));
            kwh >>= 'A';
            for (size_t j = k - 1; j < seq.size(); ++j) {
                kwh<<= seq[j];
                VERIFY(index.valid(kwh));
                if (index.get_value(kwh).count != -1u) {
                    index.get_raw_value_reference(kwh).count += 1;
                }
            }
        }

        return 0;
    }

    void ProcessCounts(IndexT &index) const {
        for (auto it = index.value_begin(); it != index.value_end(); ++it) {
            if (it->count > 1) {
                it->count = -1u;
            } else {
                it->count = 0;
            }
        }
    }

    template<class Streams>
    size_t FindRepeats(IndexT &index, Streams &streams) const {
        INFO("Collecting k-mer coverage information from reads, this takes a while.");
        unsigned nthreads = (unsigned) streams.size();
        streams.reset();
        for (size_t i = 0; i < nthreads; ++i) {
            FillCoverageFromStream(streams[i], index);
            ProcessCounts(index);
        }

        return 0;
    }

 public:

    template<class Streams>
    size_t BuildIndexFromStream(IndexT &index,
                                Streams &streams,
                                io::SingleStream* contigs_stream = 0) const {
        base::BuildIndexFromStream(index, streams, contigs_stream);

        return FindRepeats(index, streams);
    }

};

template<class Index>
struct CountIndexHelper {
    typedef Index IndexT;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::KMerIdx KMerIdx;
    typedef typename IndexT::traits_t traits_t;
//    typedef typename IndexT::IdType IdType;
    typedef DeBruijnStreamKMerIndexBuilder<Kmer, IndexT> DeBruijnStreamKMerIndexBuilderT;
    typedef RepeatSearchingIndexBuilder<DeBruijnStreamKMerIndexBuilderT> RepeatSearchingIndexBuilderT;
};

class RandNucl {
    unsigned seed_;
    boost::mt19937 rand_engine_;
    boost::uniform_int<> rand_dist_;
    boost::variate_generator<boost::mt19937&, boost::uniform_int<>> rand_nucl_;

public:

    RandNucl(unsigned seed) :
        seed_(seed),
        rand_engine_(seed_),
        rand_dist_(0, 3),
        rand_nucl_(rand_engine_, rand_dist_) {

    }

    char operator()() {
        return nucl((char) rand_nucl_());
    }
};

class RepeatMasker : public io::SequenceModifier {
private:
    typedef RtSeq Kmer;
    typedef KeyIteratingMap<Kmer, Count, kmer_index_traits<Kmer>, SimpleStoring> KmerCountIndex;
    typedef typename KmerCountIndex::KeyWithHash KeyWithHash;
    typedef KmerCountIndex::KMerIdx KmerIdx;

    size_t k_;

    RandNucl& rand_nucl_;

    KmerCountIndex index_;
    //todo maybe remove mutable? will need removing const from Modify


    bool IsRepeat(const KeyWithHash& kwh) const {
        return index_.get_value(kwh).count == -1u;
    }

    template<class S>
    const vector<Range> RepeatIntervals(const S& s) const {
        vector<Range> answer;
        answer.push_back(Range(0, 0));
        KeyWithHash kwh = index_.ConstructKWH(Kmer(k_, s) >> 'A');
        for (size_t i = k_ - 1; i < s.size(); ++i) {
            kwh <<= s[i];
            if (IsRepeat(kwh)) {
                if (i <= answer.back().end_pos) {
                    answer.back().end_pos = i + k_;
                } else {
                    answer.push_back(Range(i, i + k_));
                }
            }
        }
        return answer;
    }

    void MaskRepeat(const Range& repeat, std::string& s) const {
        TRACE("Masking repeat of size " << repeat.size() << " " << repeat);
        TRACE("Old sequence " << s.substr(repeat.start_pos, repeat.size()));
        for (size_t i = repeat.start_pos; i < repeat.end_pos; ++i) {
            s[i] = rand_nucl_();
        }
        TRACE("New sequence " << s.substr(repeat.start_pos, repeat.size()));
    }

    void MaskRepeats(const vector<Range>& rep_int, std::string& s) const {
        TRACE("Masking " << rep_int.size() << " repeat ranges in sequence of length " << s.size());
        for (Range r : rep_int) {
            MaskRepeat(r, s);
        }
    }

public:
    RepeatMasker(size_t k, RandNucl& rand_nucl, const std::string& workdir) :
        k_(k),
        rand_nucl_(rand_nucl),
        index_(k_, workdir) {
    }

    template<class Streams>
    size_t FindRepeats(Streams streams) {
        INFO("Looking for repetitive " << k_ << "-mers");
        CountIndexHelper<KmerCountIndex>::RepeatSearchingIndexBuilderT().BuildIndexFromStream(index_, streams);
        size_t rep_kmer_cnt = 0;
        for (auto it = index_.value_cbegin(); it != index_.value_cend(); ++it) {
            if (it->count == -1u) {
                rep_kmer_cnt++;
            } else {
                VERIFY(it->count == 0);
            }
        }
        INFO("Found " << rep_kmer_cnt << " repetitive " << k_ << "-mers");
        return rep_kmer_cnt;
    }

    /*virtual*/Sequence Modify(const Sequence& s) {
        if (s.size() < k_)
            return s;
        string str = s.str();
        MaskRepeats(RepeatIntervals(s), str);
        return Sequence(str);
    }

private:
    DECL_LOGGER("RepeatMasker")
};

template<class Stream1, class Stream2>
void Transfer(Stream1& s1, Stream2& s2) {
    typename Stream1::ReadT r;
    while (!s1.eof()) {
        s1 >> r;
        s2 << r;
    }
}

inline ContigStreams OpenStreams(const string& root,
        const vector<string>& filenames, bool add_rc) {
    ContigStreams streams;
    for (auto filename : filenames) {
        DEBUG("Opening stream from " << root << filename);
        ContigStreamPtr reader = make_shared<io::FileReadStream>(root + filename);
        if (add_rc)
            reader = io::RCWrap<Contig>(reader);
        streams.push_back(reader);
    }
    return streams;
}

inline void SaveStreams(ContigStreams streams, const vector<string>& suffixes,
        const string& out_root) {
    make_dir(out_root);

    streams.reset();
    for (size_t i = 0; i < streams.size(); ++i) {
        VERIFY(!suffixes[i].empty());
        string output_filename = out_root + suffixes[i];
        io::osequencestream ostream(output_filename);
        Transfer(streams[i], ostream);
    }
}

inline void ModifyAndSave(shared_ptr<io::SequenceModifier> modifier, ContigStreams streams, const vector<string>& suffixes,
        const string& out_root) {
    ContigStreams modified;
    for (size_t i = 0; i < streams.size(); ++i) {
        modified.push_back(make_shared<io::ModifyingWrapper<Contig>>(streams.ptr_at(i), modifier));
    }
    SaveStreams(modified, suffixes, out_root);
}

inline bool MaskRepeatsIteration(size_t k, const string& input_dir, const vector<string>& suffixes, const string& output_dir, RandNucl& rand_nucl) {
    shared_ptr<RepeatMasker> masker_ptr = make_shared<RepeatMasker>(k, rand_nucl, "tmp");
    INFO("Opening streams in " << input_dir);
    bool repeats_found = masker_ptr->FindRepeats(OpenStreams(input_dir, suffixes, true));
    if (repeats_found) {
        INFO("Repeats found");
        INFO("Modifying and saving streams to " << output_dir);
        ModifyAndSave(masker_ptr, OpenStreams(input_dir, suffixes, false), suffixes, output_dir);
    } else {
        INFO("No repeats found");
    }
    return !repeats_found;
}

//inline bool MaskRepeats(const string& input_dir, const vector<string>& suffixes, size_t max_iter_count, const string& work_dir) {
//    size_t iter = 0;
//    bool no_repeats = false;
//    while (iter <= max_iter_count) {
//        string out_dir = input_dir + std::to_string(iter) + "/";
//        make_dir(out_dir);
//        no_repeats = MaskRepeatsIteration(input_dir, suffixes, out_dir);
//        if (no_repeats) {
//            break;
//        }
//        ++iter;
//    }
//    if (no_repeats) {
//        string out_dir = input_dir + "masked/";
//        make_dir(out_dir);
//        ModifyAndSave(make_shared<io::TrivialModifier>(),
//                      OpenStreams(input_dir + "/" + std::to_string(iter) + "/", suffixes,
//                                  out_dir));
//    } else {
//        WARN("Failed to mask repeats in " << max_iter_count << " iterations");
//    }
//    return no_repeats;ContigStreamsPtr
//}

inline bool MaskRepeats(size_t k, ContigStreams input_streams, const vector<string>& suffixes, size_t max_iter_count, const string& work_dir) {
    size_t iter = 0;
    RandNucl rand_nucl(239);
    bool no_repeats = false;
    string input_dir = work_dir + "init/";
    make_dir(input_dir);
    SaveStreams(input_streams, suffixes, input_dir);
    while (iter <= max_iter_count) {
        INFO("------------------------");
        INFO("Iteration " << iter);
        string out_dir = work_dir + std::to_string(iter) + "/";
        make_dir(out_dir);
        no_repeats = MaskRepeatsIteration(k, input_dir, suffixes, out_dir, rand_nucl);
        if (no_repeats) {
            INFO("No repeats found");
            break;
        }
        input_dir = out_dir;
        ++iter;
    }
    if (no_repeats) {
        INFO("Repeats succesfully masked in " << iter << " iterations");
        string out_dir = work_dir + "masked/";
        make_dir(out_dir);
        SaveStreams(OpenStreams(input_dir, suffixes, false),
                    suffixes, out_dir);
    } else {
        WARN("Failed to mask repeats in " << max_iter_count << " iterations");
    }
    return no_repeats;
}

//void FillBagForStrand(const Sequence& strand,
//        map<Seq<k>>& bag) {
//    if (strand.size() < k)
//        return;
//    Seq<k> kmer(strand);
//    kmer >> 'A';
//    for (size_t i = k - 1; i < strand.size(); ++i) {
//        kmer = kmer << strand[i];
//        bag[kmer] += 1;
//    }
//}
//
//template<size_t k>
//void FillRepeats(const Sequence& genome,
//        set<Seq<k>, typename Seq<k>::less2>& repeats) {
//    map<Seq<k>, size_t, typename Seq<k>::less2> bag;
//
//    FillBagForStrand(genome, bag);
//    FillBagForStrand(!genome, bag);
//
//    for (auto it = bag.begin(); it != bag.end(); ++it) {
//        if (it->second > 1)
//            repeats.insert(it->first);
//    }
//}
//
//template<size_t k>
//void FillRepeats(const vector<Sequence>& assembly,
//        set<Seq<k>, typename Seq<k>::less2>& repeats) {
//    map<Seq<k>, size_t, typename Seq<k>::less2> bag;
//
//    for (auto it = assembly.begin(); it != assembly.end(); ++it) {
//        FillBagForStrand(*it, bag);
//        FillBagForStrand(!(*it), bag);
//    }
//
//    for (auto it = bag.begin(); it != bag.end(); ++it) {
//        if (it->second > 1)
//            repeats.insert(it->first);
//    }
//}

//template<size_t k>
//class RepeatCleaner {
//    typedef Seq<k> Kmer;
//    typedef set<Kmer, typename Kmer::less2> Repeats;
//
//    void MarkPositions(size_t start, size_t end, vector<bool>& answer) {
//        for (size_t i = start; i < end; ++i) {
//            answer[i] = true;
//        }
//    }
//
//    void MarkRepeatNucls(const Sequence& s, const Repeats& repeats, vector<bool>& answer) {
////        vector<bool> answer(s.size(), false);
//        Kmer kmer(s);
//        kmer = kmer >> 'A';
//        for (size_t i = k - 1 ; i < s.size(); ++i) {
//            kmer = kmer << s[i];
//            if (repeats.count(kmer) > 0) {
//                MarkPositions(i - k + 1, i + 1, answer);
//            }
//        }
//    }
//
//    void MarkShortIslands(, size_t threshold = k) {
//
//    }
//
//public:
//
//};

//template<size_t k>
//Sequence ClearGenome(const Sequence& genome,
//        const set<Seq<k>, typename Seq<k>::less2>& repeats) {
//    INFO("Clearing genome");
//    if (genome.size() < k)
//        return genome;
//
//    string answer;
//    for (size_t i = 0; i < k - 1; ++i) {
//        answer += nucl(genome[i]);
//    }
//    //intervals of kmers that should be resolved afterwards
//    vector<Range> repeat_intervals;
//    Seq<k> kmer(genome);
//    size_t curr_pos = 0;
//    //curr_pos + k - next nucl pos
//    bool changed = false;
//    while (curr_pos + k != genome.size()) {
//        size_t int_start = curr_pos;
//        while (repeats.count(kmer) > 0 && curr_pos + k < genome.size()) {
//            kmer = kmer << genome[curr_pos + k];
//            curr_pos++;
//            changed = true;
//        }
//        if (int_start != curr_pos)
//            repeat_intervals.push_back(Range(int_start, curr_pos));
//
//        if (curr_pos + k == genome.size())
//            break;
//
//        while (repeats.count(kmer) == 0 && curr_pos + k < genome.size()) {
//            answer += nucl(kmer[k - 1]);
//            kmer = kmer << genome[curr_pos + k];
//            curr_pos++;
//        }
//    }
//    if (changed) {
//        INFO("Genome was changed during cleaning");
//    } else {
//        INFO("Genome wasn't changed during cleaning");
//    }
//    return Sequence(answer);
//}
//
//template<size_t k>
//Sequence ClearGenome(const Sequence& genome) {
//    INFO("Clearing genome of repeats");
//
//    set<Seq<k>, typename Seq<k>::less2> repeats;
//    INFO("Filling set of repeats");
//    FillRepeats<k>(genome, repeats);
//    INFO("Clearing genome");
//    return ClearGenome<k>(genome, repeats);
//}
//
////todo bad strategy for assembly cleaning
//template<size_t k>
//pair<Sequence, vector<Sequence>> Clear(const Sequence& genome,
//        const vector<Sequence>& assembly) {
//    INFO("Clearing genome of repeats");
//
//    set<Seq<k>, typename Seq<k>::less2> repeats;
//    INFO("Filling set of repeats");
//    FillRepeats<k>(genome, repeats);
////    for (auto it = assembly.begin(); it != assembly.end(); ++it) {
////        FillRepeats(*it, repeats);
////    }
//    INFO("Clearing genome");
//    Sequence new_genome = ClearGenome<k>(genome, repeats);
//    INFO("Clearing assembly");
//    vector<Sequence> new_assembly;
//    for (auto it = assembly.begin(); it != assembly.end(); ++it) {
//        new_assembly.push_back(ClearGenome<k>(*it, repeats));
//    }
//    return make_pair(new_genome, new_assembly);
//}

//template<size_t k>
//void Clear(const string& genome_in, const string& genome_out
//        , const string& assembly_in, const string& assembly_out) {
//    INFO("Clearing genome of repeats");
//    pair<Sequence, vector<Sequence>> cleared
//        = Clear<k>(ReadGenome(genome_in), ReadContigs(assembly_in));
//    io::ofastastream genome_out
//
//}

//template<size_t k>
//void Clear(const string& in, const string& out) {
//    io::Reader in_stream(in);
//    set<Seq<k>, typename Seq<k>::less2> repeats;
//    FillRepeats<k>(AllSequences(in_stream), repeats);
//    in_stream.reset();
//    io::osequencestream out_stream(out);
//    io::SingleRead contig;
//    while (!in_stream.eof()) {
//        in_stream >> contig;
//        Sequence cleared = ClearGenome<k>(contig.sequence(), repeats);
//        out_stream << io::SingleRead(contig.name(), cleared.str());
//    }
//}
//
//template<size_t k>
//pair<Sequence, Sequence> ClearGenomes(const pair<Sequence, Sequence>& genomes) {
//    INFO("Clearing genomes from repeats");
//
//    set<Seq<k>, typename Seq<k>::less2> repeats;
//    INFO("Filling set of repeats");
//    FillRepeats<k>(genomes.first, repeats);
//    FillRepeats<k>(genomes.second, repeats);
//    INFO("Clearing genomes");
//    return make_pair(ClearGenome<k>(genomes.first, repeats),
//            ClearGenome<k>(genomes.second, repeats));
//}
//
//template<size_t k>
//pair<Sequence, Sequence> TotallyClearGenomes(
//        const pair<Sequence, Sequence>& genomes) {
//    static const size_t iter_count = 1;
//    pair<Sequence, Sequence> tmp = genomes;
//    for (size_t i = 0; i < iter_count; ++i) {
//        INFO("Cleaning iteration " << i);
//        tmp = ClearGenomes<k>(tmp);
//    }
//    return tmp;
//}
//
//template<size_t k>
//bool CheckNoRepeats(const Sequence& genome) {
//    set<Seq<k>, typename Seq<k>::less2> repeats;
//    FillRepeats<k>(genome, repeats);
//    return repeats.empty();
//}



}
