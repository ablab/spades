template<class Data>
using Library = io::SequencingLibrary<Data>;

namespace llvm { namespace yaml {
template <class Data>
struct SequenceTraits<std::vector<Library<Data> >>  {
    static size_t size(IO &, std::vector<Library<Data> > &seq) {
        return seq.size();
    }
    static Library<Data>&
    element(IO &, std::vector<Library<Data>> &seq, size_t index) {
        if (index >= seq.size())
            seq.resize(index+1);
        return seq[index];
    }
};

template<class Data>
void MappingTraits<Library<Data>>::mapping(yaml::IO &io, Library<Data> &lib) {
    lib.yamlize(io);
}

template<class Data>
StringRef MappingTraits<Library<Data>>::validate(yaml::IO &io, Library<Data> &lib) {
    // We use such crazy API for validate() since we don't want to pull
    // llvm::StringRef into library.hpp.
    llvm::StringRef res;
    lib.validate(io, res);

    return res;
}

template<class Data>
void MappingTraits<io::DataSet<Data>>::mapping(yaml::IO &io, io::DataSet<Data> &ds) {
    ds.yamlize(io);
}
}}

template<class Data>
void io::SequencingLibrary<Data>::yamlize(llvm::yaml::IO &io) {
    // First, load the "common stuff"
    SequencingLibraryBase::yamlize(io);
    io.mapOptional("data", data_);
}

template<class Data>
void io::SequencingLibrary<Data>::validate(llvm::yaml::IO &io, llvm::StringRef &res) {
    // Simply ask base class to validate for us
    SequencingLibraryBase::validate(io, res);
}

template<class Data>
void io::DataSet<Data>::yamlize(yaml::IO &io) {
    llvm::yaml::yamlize(io, libraries_, true);
}

template<class Data>
void io::DataSet<Data>::save(const std::string &filename) {
    std::error_code EC;
    llvm::raw_fd_ostream ofs(filename, EC, llvm::sys::fs::OpenFlags::F_Text);
    llvm::yaml::Output yout(ofs);
    yout << libraries_;
}

template<class Data>
void io::DataSet<Data>::load(const std::string &filename) {
    ErrorOr<std::unique_ptr<MemoryBuffer>> Buf = MemoryBuffer::getFile(filename);
    if (!Buf) {
        std::cerr << std::string("Failed to load file ") + filename;
        throw;
    }

    yaml::Input yin(*Buf.get());
    yin >> libraries_;

    if (yin.error()) {
        std::cerr << std::string("Failed to load file ") + filename;
        throw;
    }
    
    std::string input_dir = fs::parent_path(filename);
    if (input_dir[input_dir.length() - 1] != '/')
        input_dir += '/';

    for (auto& lib : libraries_)
        lib.update_relative_reads_filenames(input_dir);
}
