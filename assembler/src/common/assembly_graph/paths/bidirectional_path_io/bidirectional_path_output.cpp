//
// Created by andrey on 20.01.17.
//

#include "bidirectional_path_output.hpp"

namespace path_extend {


string path_extend::ContigWriter::ToFASTGPathFormat(const BidirectionalPath &path) const {
    if (path.Empty())
        return "";
    string res = ids_.at(path.Front()).short_id_;
    for (size_t i = 1; i < path.Size(); ++i) {
        if (g_.EdgeEnd(path[i - 1]) != g_.EdgeStart(path[i])) {
            res += ";\n" + ids_.at(path[i]).short_id_;
        } else {
            res += "," + ids_.at(path[i]).short_id_;
        }
    }
    return res;
}

void path_extend::ContigWriter::WriteScaffolds(const vector<ScaffoldInfo> &scaffold_storage, const string &fn) const {
    io::osequencestream_simple oss(fn);
    std::ofstream os_fastg;

    for (const auto& scaffold_info : scaffold_storage) {
        oss.set_header(scaffold_info.name);
        oss << scaffold_info.sequence;
    }
}

void path_extend::ContigWriter::WritePathsFastg(const vector<ScaffoldInfo> &scaffold_storage, const string &fn) const {
    std::ofstream os(fn);
    for (const auto& scaffold_info : scaffold_storage) {
        os << scaffold_info.name << endl;
        os << ToFASTGPathFormat(*scaffold_info.path) << endl;
        os << scaffold_info.name << "'" << endl;
        os << ToFASTGPathFormat(*scaffold_info.path->GetConjPath()) << endl;
    }
}

void path_extend::ContigWriter::OutputPaths(const PathContainer &paths,
                                            const string &filename_base,
                                            bool write_fastg) const {
    vector<ScaffoldInfo> storage;

    ScaffoldSequenceMaker scaffold_maker(g_);
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        BidirectionalPath* path = iter.get();
        if (path->Length() <= 0)
            continue;
        string path_string = scaffold_maker.MakeSequence(*path);
        if (path_string.length() >= g_.k()) {
            storage.emplace_back(path_string, path);
        }
    }

    //sorting by length and coverage
    std::sort(storage.begin(), storage.end(), [] (const ScaffoldInfo &a, const ScaffoldInfo &b) {
        if (a.length() == b.length())
            return math::gr(a.coverage(), b.coverage());
        return a.length() > b.length();
    });

    name_generator_->Preprocess(paths);
    for (size_t i = 0; i < storage.size(); ++i) {
        storage[i].name = name_generator_->MakeContigName(i+1, storage[i]);
    }

    INFO("Writing contigs to " << filename_base);
    if (write_fastg)
        WritePathsFastg(storage, filename_base + ".paths");

    WriteScaffolds(storage, filename_base + ".fasta");

    DEBUG("Contigs written");
}

}