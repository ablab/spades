//
// Created by andrey on 20.01.17.
//

#include "bidirectional_path_output.hpp"

namespace path_extend {

void path_extend::ContigWriter::OutputPaths(const PathContainer &paths, const vector<PathsWriterT>& writers) const {
    ScaffoldStorage storage;

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

    for (auto& writer : writers) {
        writer(storage);
    }
    DEBUG("Contigs written");
}

}