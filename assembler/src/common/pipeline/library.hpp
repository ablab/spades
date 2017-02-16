//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __IO_LIBRARY_HPP__
#define __IO_LIBRARY_HPP__

#include "common/adt/chained_iterator.hpp"
#include "common/adt/iterator_range.hpp"

#include <boost/iterator/iterator_facade.hpp>

#include <string>
#include <vector>

// Forward decls for YAML API
namespace llvm { namespace yaml { class IO; template<typename T> struct MappingTraits; } }
namespace llvm { class StringRef;  }

namespace io {

enum class LibraryType {
    SingleReads,
    SangerReads,
    PacBioReads,
    NanoporeReads,
    PairedEnd,
    HQMatePairs,
    MatePairs,
    TrustedContigs,
    TSLReads,
    PathExtendContigs,
    UntrustedContigs

};

static std::vector<LibraryType> LibraryPriotity = {
    LibraryType::SingleReads,
    LibraryType::SangerReads,
    LibraryType::PacBioReads,
    LibraryType::NanoporeReads,
    LibraryType::PairedEnd,
    LibraryType::HQMatePairs,
    LibraryType::MatePairs,
    LibraryType::TrustedContigs,
    LibraryType::TSLReads,
    LibraryType::PathExtendContigs,
    LibraryType::UntrustedContigs
};

enum class LibraryOrientation {
  FR,
  FF,
  RF,
  RR,
  Undefined
};

class SequencingLibraryBase {
public:
    class paired_reads_iterator :
            public boost::iterator_facade<paired_reads_iterator,
                                          std::pair<std::string, std::string>,
                                          boost::forward_traversal_tag,
                                          std::pair<std::string, std::string> > {

        typedef std::vector<std::string>::const_iterator inner_iterator;

    public:
        paired_reads_iterator(inner_iterator left, inner_iterator right)
                : left_(left), right_(right){}

    private:
        friend class boost::iterator_core_access;

        void increment() { ++left_; ++right_; }
        bool equal(const paired_reads_iterator &other) const {
            return this->left_ == other.left_ && this->right_ == other.right_;
        }
        std::pair<std::string, std::string> dereference() const {
            return std::make_pair(*left_, *right_);
        }

        inner_iterator left_;
        inner_iterator right_;
    };

    typedef chained_iterator<std::vector<std::string>::const_iterator> single_reads_iterator;

    SequencingLibraryBase()
            : type_(LibraryType::PairedEnd), orientation_(LibraryOrientation::FR) {}

    // YAML API. Public because we cannot have template friend class.
    void yamlize(llvm::yaml::IO &io);
    void validate(llvm::yaml::IO &io, llvm::StringRef &res);

    LibraryType type() const { return type_; }
    void set_type(LibraryType type) { type_ = type; }
    LibraryOrientation orientation() const { return orientation_; }
    void set_orientation(LibraryOrientation orientation) { orientation_ = orientation; }

    void clear() {
        left_paired_reads_.clear();
        right_paired_reads_.clear();
        single_reads_.clear();
    }

    void update_relative_reads_filenames(const std::string &input_dir);

    void push_back_single(const std::string &reads) {
        single_reads_.push_back(reads);
    }

    void push_back_paired(const std::string &left, const std::string &right) {
        left_paired_reads_.push_back(left);
        right_paired_reads_.push_back(right);
    }

    paired_reads_iterator paired_begin() const {
        return paired_reads_iterator(left_paired_reads_.begin(), right_paired_reads_.begin());
    }
    paired_reads_iterator paired_end() const {
        return paired_reads_iterator(left_paired_reads_.end(), right_paired_reads_.end());
    }

    adt::iterator_range<paired_reads_iterator> paired_reads() const {
        return adt::make_range(paired_begin(), paired_end());
    }

    single_reads_iterator reads_begin() const {
        // NOTE: We have a contract with single_end here. Single reads always go last!
        single_reads_iterator res(left_paired_reads_.begin(), left_paired_reads_.end());
        res.join(right_paired_reads_.begin(), right_paired_reads_.end());
        res.join(single_reads_.begin(), single_reads_.end());
        
        return res;
    }
    single_reads_iterator reads_end() const {
        // NOTE: Do not forget about the contract with single_begin here!
        return single_reads_iterator(single_reads_.end(), single_reads_.end());
    }

    adt::iterator_range<single_reads_iterator> reads() const {
        return adt::make_range(reads_begin(), reads_end());
    }

    single_reads_iterator single_begin() const {
        return single_reads_iterator(single_reads_.begin(), single_reads_.end());
    }
    single_reads_iterator single_end() const {
        // NOTE: Do not forget about the contract with single_begin here!
        return single_reads_iterator(single_reads_.end(), single_reads_.end());
    }

    adt::iterator_range<single_reads_iterator> single_reads() const {
        return adt::make_range(single_begin(), single_end());
    }

    bool is_graph_contructable() const {
        return type_ == io::LibraryType::PairedEnd ||
               type_ == io::LibraryType::SingleReads ||
               type_ == io::LibraryType::HQMatePairs;
    }

    bool is_bwa_alignable() const {
        return type_ == io::LibraryType::MatePairs;
    }

    bool is_mismatch_correctable() const {
        return is_graph_contructable();
    }

//    bool is_binary_covertable() {
//        return is_graph_contructable() || is_mismatch_correctable() || is_paired();
//    }

    bool is_paired() const {
        return type_ == io::LibraryType::PairedEnd ||
               type_ == io::LibraryType::MatePairs ||
               type_ == io::LibraryType::HQMatePairs;
    }

    bool is_mate_pair() const {
        return type_ == io::LibraryType::MatePairs ||
               type_ == io::LibraryType::HQMatePairs;
    }

    static bool is_contig_lib(LibraryType type) {
        return type == io::LibraryType::TrustedContigs ||
               type == io::LibraryType::UntrustedContigs ||
               type == io::LibraryType::PathExtendContigs;
    }

    static bool is_long_read_lib(LibraryType type) {
        return type == io::LibraryType::PacBioReads ||
               type == io::LibraryType::SangerReads ||
               type == io::LibraryType::NanoporeReads ||
               type == io::LibraryType::TSLReads;
    }

    bool is_contig_lib() const {
        return is_contig_lib(type_);
    }

    bool is_long_read_lib() const {
        return is_long_read_lib(type_);
    }

    bool is_repeat_resolvable() const {
        return is_paired() ||
               is_long_read_lib() ||
               is_contig_lib();
    }

    //hybrid libraries are used to close gaps in the graph during their alignment
    bool is_hybrid_lib() const {
        return is_long_read_lib() ||
               //comment next line to switch alignment method for trusted contigs
               type_ == io::LibraryType::TrustedContigs ||
               type_ == io::LibraryType::UntrustedContigs;
    }

private:
    LibraryType type_;
    LibraryOrientation orientation_;

    std::vector<std::string> left_paired_reads_;
    std::vector<std::string> right_paired_reads_;
    std::vector<std::string> single_reads_;
};

struct NoData {};

template<class Data = NoData>
class SequencingLibrary: public SequencingLibraryBase {
public:
    const Data& data() const {
        return data_;
    }
    Data& data() {
        return data_;
    }

    void yamlize(llvm::yaml::IO &io);
    void validate(llvm::yaml::IO &io, llvm::StringRef &res);

private:
    Data data_;
};

// Just convenient wrapper to "unwrap" the iterators over libraries.
template<class Data = NoData>
class DataSet {
public:
    typedef SequencingLibrary<Data> Library;
    typedef std::vector<Library> LibraryStorage;

public:
    typedef typename LibraryStorage::iterator iterator;
    typedef typename LibraryStorage::const_iterator const_iterator;
    typedef chained_iterator<typename Library::single_reads_iterator> single_reads_iterator;
    typedef chained_iterator<typename Library::paired_reads_iterator> paired_reads_iterator;

    DataSet() {}
    explicit DataSet(const std::string &path) { load(path); }

    void load(const std::string &filename);
    void save(const std::string &filename);

    void clear() { libraries_.clear(); }
    void push_back(const Library &lib) {
        libraries_.push_back(lib);
    }
    Library& operator[](size_t n) { return libraries_[n]; }
    const Library& operator[](size_t n) const { return libraries_[n]; }
    size_t lib_count() const { return libraries_.size(); }

    iterator library_begin() { return libraries_.begin(); }
    const_iterator library_begin() const { return libraries_.begin(); }
    iterator begin() { return libraries_.begin(); }
    const_iterator begin() const { return libraries_.begin(); }

    iterator library_end() { return libraries_.end(); }
    const_iterator library_end() const { return libraries_.end(); }
    iterator end() { return libraries_.end(); }
    const_iterator end() const { return libraries_.end(); }

    adt::iterator_range<iterator> libraries() {
        return adt::make_range(library_begin(), library_end());
    }
    adt::iterator_range<const_iterator> libraries() const {
        return adt::make_range(library_begin(), library_end());
    }

    single_reads_iterator reads_begin() const {
        auto it = libraries_.begin();
        single_reads_iterator res(it->reads_begin(), it->reads_end());
        ++it;
        for (auto end = libraries_.end(); it != end; ++it)
            res.join(it->reads_begin(), it->reads_end());

        return res;
    }
    single_reads_iterator reads_end() const {
        return single_reads_iterator(libraries_.back().reads_end(), libraries_.back().reads_end());
    }
    adt::iterator_range<single_reads_iterator> reads() const {
        return adt::make_range(reads_begin(), reads_end());
    }

    single_reads_iterator single_begin() const {
        auto it = libraries_.begin();
        single_reads_iterator res(it->single_begin(), it->single_end());
        ++it;
        for (auto end = libraries_.end(); it != end; ++it)
            res.join(it->single_begin(), it->single_end());

        return res;
    }
    single_reads_iterator single_end() const {
        return single_reads_iterator(libraries_.back().single_end(), libraries_.back().single_end());
    }
    adt::iterator_range<single_reads_iterator> single_reads() const {
        return adt::make_range(single_begin(), single_end());
    }

    paired_reads_iterator paired_begin() const {
        auto it = libraries_.begin();
        paired_reads_iterator res(it->paired_begin(), it->paired_end());
        ++it;
        for (auto end = libraries_.end(); it != end; ++it)
            res.join(it->paired_begin(), it->paired_end());

        return res;
    }
    paired_reads_iterator paired_end() const {
        return paired_reads_iterator(libraries_.back().paired_end(), libraries_.back().paired_end());
    }

    adt::iterator_range<paired_reads_iterator> paired_reads() const {
        return adt::make_range(paired_begin(), paired_end());
    }

private:
    LibraryStorage libraries_;
};

}

namespace llvm { namespace yaml {
template <>
struct MappingTraits<io::SequencingLibraryBase> {
    static void mapping(llvm::yaml::IO &io, io::SequencingLibraryBase &lib);
    static StringRef validate(llvm::yaml::IO &io, io::SequencingLibraryBase &lib);
};

template <class Data>
struct MappingTraits<io::SequencingLibrary<Data> > {
    static void mapping(llvm::yaml::IO &io, io::SequencingLibrary<Data> &lib);
    static StringRef validate(llvm::yaml::IO &io, io::SequencingLibrary<Data> &lib);
};

}}

#endif // __IO_LIBRARY_HPP__
