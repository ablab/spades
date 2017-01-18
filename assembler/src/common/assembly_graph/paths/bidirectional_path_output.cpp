//
// Created by andrey on 20.01.17.
//

#include "bidirectional_path_output.hpp"
#include "modules/path_extend/pe_utils.hpp"
#include "io/reads/osequencestream.hpp"


void path_extend::TranscriptToGeneJoiner::MakeSet(size_t x) {
    parents_[x] = x;
    ranks_[x] = 0;
}

void path_extend::TranscriptToGeneJoiner::JoinTrees(size_t x, size_t y) {
    x = FindTree(x);
    y = FindTree(y);
    if (x != y) {
        if (ranks_[x] < ranks_[y])
            parents_[x] = y;
        else
            parents_[y] = x;
        if (ranks_[x] == ranks_[y])
            ++ranks_[x];
    }
}

void path_extend::TranscriptToGeneJoiner::Init(const PathContainer &paths) {
    DEBUG("Initializing parents and ranks");
    parents_.resize(paths.size());
    ranks_.resize(paths.size());

    size_t path_num = 0;
    for (auto iter = paths.begin(); iter != paths.end(); ++iter, ++path_num) {
        path_id_[iter.get()] = path_num;
        path_id_[iter.getConjugate()] = path_num;
        MakeSet(path_num);
    }

    DEBUG("Initialized parents and ranks");

    VERIFY_MSG(path_num == paths.size(), "Path Num " << path_num << " Size " << paths.size())
}

size_t path_extend::TranscriptToGeneJoiner::FindTree(size_t x) {
    size_t parent;
    if (x == parents_[x]) {
        parent = x;
    }
    else {
        parents_[x] = FindTree(parents_[x]);
        parent = parents_[x];
    }
    return parent;
}

size_t path_extend::TranscriptToGeneJoiner::GetPathId(BidirectionalPath *path) {
    return path_id_[path];
}

void path_extend::TranscriptToGeneJoiner::Construct(const PathContainer &paths) {
    Init(paths);

    GraphCoverageMap edges_coverage(g_, paths);

    DEBUG("Union trees");
    //For all edges in coverage map
    for (auto iterator = edges_coverage.begin(); iterator != edges_coverage.end(); ++iterator) {
        //Select a path covering an edge
        EdgeId edge = iterator->first;
        GraphCoverageMap::MapDataT *edge_paths = iterator->second;

        if (g_.length(edge) > min_edge_len_ && edge_paths->size() > 1) {
            DEBUG("Long edge " << edge.int_id() << " Paths " << edge_paths->size());
            //For all other paths covering this edge join then into single gene with the first path
            for (auto it_edge = ++edge_paths->begin(); it_edge != edge_paths->end(); ++it_edge) {
                size_t first = path_id_[*edge_paths->begin()];
                size_t next = path_id_[*it_edge];
                DEBUG("Edge " << edge.int_id() << " First " << first << " Next " << next);

                JoinTrees(first, next);
            }
        }
    }
}

string path_extend::ContigWriter::ToString(const BidirectionalPath& path) const {
    stringstream ss;
    if (path.IsInterstrandBulge() && path.Size() == 1) {
        ss << constructor_.construct(path.Back()).first.substr(k_, g_.length(path.Back()) - k_);
        return ss.str();
    }

    if (!path.Empty()) {
        ss << constructor_.construct(path[0]).first.substr(0, k_);
    }


    size_t i = 0;
    while (i < path.Size()) {
        int gap = i == 0 ? 0 : path.GapAt(i);
        if (gap > (int) k_) {
            for (size_t j = 0; j < gap - k_; ++j) {
                ss << "N";
            }
            ss << constructor_.construct(path[i]).first;
        }
        else {
            int overlapLen = (int) k_ - gap;
            if (overlapLen >= (int) g_.length(path[i]) + (int) k_) {
                overlapLen -= (int) g_.length(path[i]) + (int) k_;
                ++i;
                //skipping overlapping edges
                while (i < path.Size() && overlapLen >= (int) g_.length(path[i]) + path.GapAt(i)) {
                    overlapLen -= (int) g_.length(path[i]) + path.GapAt(i);
                    ++i;
                }
                if (i == path.Size()) {
                    break;
                }

                overlapLen = overlapLen + (int) k_ - path.GapAt(i);
                if(overlapLen < 0) {
                    for (int j = 0; j < abs(overlapLen); ++j) {
                        ss << "N";
                    }
                    overlapLen = 0;
                }
            }
            auto temp_str = g_.EdgeNucls(path[i]).Subseq(overlapLen).str();
            if(i != path.Size() - 1) {
                for(size_t j = 0 ; j < path.TrashPreviousAt(i + 1); ++j) {
                    temp_str.pop_back();
                    if(temp_str.size() == 0) {
                        break;
                    }
                }
            }
            ss << temp_str;
        }
        ++i;
    }
    return ss.str();
}

string path_extend::ContigWriter::ToFASTGString(const BidirectionalPath& path) const {
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

void path_extend::ContigWriter::WritePathsToFASTA(const PathContainer &paths,
                                                  const string &filename_base,
                                                  bool write_fastg, size_t long_edge_threshold) const {

    vector<IOContigStorage> storage;
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        BidirectionalPath* path = iter.get();
        if (path->Length() <= 0)
            continue;
        string path_string = ToString(*path);
        if (path_string.length() >= g_.k()) {
            storage.emplace_back(path_string, path);
        }
    }
    std::sort(storage.begin(), storage.end(), IOContigStorageGreater());

    INFO("Writing contigs to " << filename_base);
    io::osequencestream_simple oss(filename_base + ".fasta");
    std::ofstream os_fastg;
    if (write_fastg)
        os_fastg.open((filename_base + ".paths").c_str());

    TranscriptToGeneJoiner transcript_joiner(g_, long_edge_threshold);
    if (pipeline_type_ == config::pipeline_type::rna) {
        transcript_joiner.Construct(paths);
    }

    size_t gene_num = 0;
    map<size_t, size_t> isoform_num;
    map<size_t, size_t> gene_ids;
    size_t i = 0;
    for (const auto& precontig : storage) {
        ++i;
        std::string contig_id = "";

        if (pipeline_type_ ==  config::pipeline_type::plasmid) {
            EdgeId e = precontig.path_->At(0);
            size_t component = c_counter_.GetComponent(e);
            contig_id = io::MakeContigComponentId(i, precontig.sequence_.length(), precontig.path_->Coverage(), component);
        }
        else if (pipeline_type_ == config::pipeline_type::rna) {
            size_t id = transcript_joiner.GetPathId(precontig.path_);
            size_t parent_id = transcript_joiner.FindTree(id);
            DEBUG("Path " << id << " Parent " << parent_id);
            if (gene_ids.find(parent_id) == gene_ids.end()) {
                gene_ids[parent_id] = gene_num;
                isoform_num[parent_id] = 0;
                gene_num++;
            }
            contig_id = io::MakeRNAContigId(i, precontig.sequence_.length(), precontig.path_->Coverage(), gene_ids[parent_id], isoform_num[parent_id]);
            isoform_num[parent_id]++;
        } else {
            contig_id = io::MakeContigId(i, precontig.sequence_.length(), precontig.path_->Coverage());
        }
        oss.set_header(contig_id);
        oss << precontig.sequence_;

        if (write_fastg) {
            os_fastg << contig_id << endl;
            os_fastg << ToFASTGString(*precontig.path_) << endl;
            os_fastg << contig_id << "'" << endl;
            os_fastg << ToFASTGString(*precontig.path_->GetConjPath()) << endl;
        }
    }

    if (write_fastg)
        os_fastg.close();
    DEBUG("Contigs written");
}

void path_extend::ContigWriter::OutputPaths(const PathContainer& paths, const string& filename_base) const {
    WritePathsToFASTA(paths, filename_base);
}

void path_extend::PathInfoWriter::WritePaths(const PathContainer &paths, const string &filename) const {
    std::ofstream oss(filename.c_str());

    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        iter.get()->Print(oss);
    }

    oss.close();
}

void path_extend::ScaffoldBreaker::SplitPath(const BidirectionalPath& path, PathContainer &result) const {
    size_t i = 0;

    while (i < path.Size()) {
        BidirectionalPath * p = new BidirectionalPath(path.graph(), path[i]);
        ++i;

        while (i < path.Size() and path.GapAt(i) <= min_gap_) {
            p->PushBack(path[i], path.GapAt(i), path.TrashPreviousAt(i), path.TrashCurrentAt(i));
            ++i;
        }

        if (i < path.Size()) {
            DEBUG("split path " << i << " gap " << path.GapAt(i));
            p->Print();
        }

        BidirectionalPath * cp = new BidirectionalPath(p->Conjugate());
        result.AddPair(p, cp);
    }
}

void path_extend::ScaffoldBreaker::Break(const PathContainer &paths, PathContainer &result) const {
    for (auto it = paths.begin(); it != paths.end(); ++it) {
        SplitPath(*it.get(), result);
    }
    result.SortByLength();
}
