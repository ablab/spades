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
        }
        else {
            res += "," + ids_.at(path[i]).short_id_;
        }
    }
    return res;
}

void path_extend::ContigWriter::OutputPaths(const PathContainer &paths,
                                                  const string &filename_base,
                                                  bool write_fastg) const {
    name_generator_->Preprocess(paths);
    IOContigStorage storage(g_, constructor_, paths);

    INFO("Writing contigs to " << filename_base);
    io::osequencestream_simple oss(filename_base + ".fasta");
    std::ofstream os_fastg;
    if (write_fastg)
        os_fastg.open((filename_base + ".paths").c_str());

    size_t i = 0;
    for (const auto& precontig : storage.Storage()) {
        ++i;
        std::string contig_id = name_generator_->MakeContigName(i, precontig);
        oss.set_header(contig_id);
        oss << precontig.sequence_;

        if (write_fastg) {
            os_fastg << contig_id << endl;
            os_fastg << ToFASTGPathFormat(*precontig.path_) << endl;
            os_fastg << contig_id << "'" << endl;
            os_fastg << ToFASTGPathFormat(*precontig.path_->GetConjPath()) << endl;
        }
    }

    if (write_fastg)
        os_fastg.close();
    DEBUG("Contigs written");
}


void path_extend::PathInfoWriter::WritePaths(const PathContainer &paths, const string &filename) const {
    std::ofstream oss(filename.c_str());

    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        iter.get()->Print(oss);
    }

    oss.close();
}

}