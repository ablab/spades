//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef PE_IO_HPP_
#define PE_IO_HPP_


#include "modules/genome_consistance_checker.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/graph_support/contig_output.hpp"
#include "assembly_graph/components/connected_component.hpp"
#include "io/reads/osequencestream.hpp"
#include <stack>
#include <algorithm>

namespace path_extend {

using namespace debruijn_graph;

struct IOContigStorage {
    std::string sequence_;
    BidirectionalPath* path_;

    IOContigStorage(const std::string& sequence, BidirectionalPath* path) :
    sequence_(sequence), path_(path) { }
};

struct IOContigStorageGreater
{
    bool operator()(const IOContigStorage &a, const IOContigStorage &b) const {
        if (a.sequence_.length() ==  b.sequence_.length())
            return math::gr(a.path_->Coverage(), b.path_->Coverage());
        return a.sequence_.length() > b.sequence_.length();
    }
};

//Finds common long edges in paths and joins them into
class TranscriptToGeneJoiner {
private:
    const Graph &g_;
    size_t min_edge_len_; //minimal length for joining transcripts into a gene

    map<BidirectionalPath *, size_t, PathComparator> path_id_; //path ids
    std::vector<size_t> parents_; //node parents in
    std::vector<size_t> ranks_; //tree depth


    void MakeSet(size_t x) {
        parents_[x] = x;
        ranks_[x] = 0;
    }

    void UnionTrees(size_t x, size_t y) {
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

    void Init(const PathContainer &paths) {
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
public:
    TranscriptToGeneJoiner(const Graph &g, size_t min_edge_len): g_(g), min_edge_len_(min_edge_len) {}

    size_t FindTree(size_t x) {
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

    size_t GetPathId(BidirectionalPath *path) {
        return path_id_[path];
    }

    void Construct(const PathContainer &paths) {
        Init(paths);

        GraphCoverageMap edges_coverage(g_, paths);

        DEBUG("Union trees");
        for (auto iterator = edges_coverage.begin(); iterator != edges_coverage.end(); ++iterator) {
            EdgeId edge = iterator->first;
            GraphCoverageMap::MapDataT *edge_paths = iterator->second;

            if (g_.length(edge) > min_edge_len_ && edge_paths->size() > 1) {
                DEBUG("Long edge " << edge.int_id() << " Paths " << edge_paths->size());
                for (auto it_edge = ++edge_paths->begin(); it_edge != edge_paths->end(); ++it_edge) {
                    size_t first = path_id_[*edge_paths->begin()];
                    size_t next = path_id_[*it_edge];
                    DEBUG("Edge " << edge.int_id() << " First " << first << " Next " << next);

                    UnionTrees(first, next);
                }
            }
        }
    }
};

class ContigWriter {
protected:
    DECL_LOGGER("PathExtendIO")

protected:
    const Graph& g_;
    ContigConstructor<Graph> &constructor_;
    size_t k_;
    map<EdgeId, ExtendedContigIdT> ids_;
    const ConnectedComponentCounter &c_counter_;
    config::pipeline_type pipeline_type_;

    //TODO: add constructor
    string ToString(const BidirectionalPath& path) const {
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

    string ToFASTGString(const BidirectionalPath& path) const {
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


public:
    ContigWriter(const Graph& g,
                 ContigConstructor<Graph> &constructor,
                 const ConnectedComponentCounter &c_counter,
                 config::pipeline_type pipeline_type) :
            g_(g), constructor_(constructor), k_(g.k()),
            ids_(), c_counter_(c_counter),
            pipeline_type_(pipeline_type) {
        MakeContigIdMap(g_, ids_, c_counter, "NODE");
    }


    void WriteEdges(const string &filename) const {
        INFO("Outputting edges to " << filename);
        io::osequencestream_with_id oss(filename);

        set<EdgeId> included;
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (included.count(*iter) == 0) {
                oss.setCoverage(g_.coverage(*iter));
                oss.setID((int) g_.int_id(*iter));
                oss << g_.EdgeNucls(*iter);

                included.insert(*iter);
                included.insert(g_.conjugate(*iter));
            }
        }
        DEBUG("Contigs written");
    }


    void WritePaths(const PathContainer &paths, const string &filename) const {
        INFO("Outputting path data to " << filename);
        std::ofstream oss;
        oss.open(filename.c_str());
        int i = 1;
        oss << paths.size() << endl;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            //oss << i << endl;
            i++;
            BidirectionalPath* path = iter.get();
            if (path->GetId() % 2 != 0) {
                path = path->GetConjPath();
            }
            oss << "PATH " << path->GetId() << " " << path->Size() << " " << path->Length() + k_ << endl;
            for (size_t j = 0; j < path->Size(); ++j) {
                oss << g_.int_id(path->At(j)) << " " << g_.length(path->At(j)) <<  " " << path->GapAt(j) <<  " " << path->TrashPreviousAt(j) <<  " " << path->TrashCurrentAt(j) << endl;
            }
            //oss << endl;
        }
        oss.close();
        DEBUG("Edges written");
    }

    void LoadPaths(PathContainer &paths, GraphCoverageMap &cover_map, const string &filename) const {
        paths.clear();
        map<size_t, EdgeId> int_ids;
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            int_ids.insert(make_pair(g_.int_id(*iter), *iter));
        }

        std::ifstream iss;
        iss.open(filename);
        size_t psize;
        iss >> psize;
        for(size_t i = 0; i < psize && !iss.eof(); ++i) {
            string s;
            size_t id;
            size_t size;
            size_t len;
            iss >> s >> id >> size >> len;
            VERIFY(s == "PATH");

            BidirectionalPath * path = new BidirectionalPath(g_);
            BidirectionalPath * conjugatePath = new BidirectionalPath(g_);
            paths.AddPair(path, conjugatePath);
            path->Subscribe(&cover_map);
            conjugatePath->Subscribe(&cover_map);
            for (size_t j = 0; !iss.eof() && j < size; ++j) {
                size_t eid;
                size_t elen;
                int gap;
                uint32_t trash_prev;
                uint32_t trash_current;

                iss >> eid >> elen >> gap >> trash_prev >> trash_current;
                Gap gap_struct(gap, trash_prev, trash_current);
                EdgeId edge = int_ids[eid];
                conjugatePath->PushBack(edge, gap_struct);
                VERIFY(g_.length(edge) == elen);
            }
            VERIFY(path->Length() + k_ == len);
        }
        VERIFY(psize == paths.size());
        iss.close();
    }

    void WritePathsToFASTA(const PathContainer &paths,
                           const string &filename_base,
                           bool write_fastg = true, size_t long_length = 300) const {

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

        TranscriptToGeneJoiner transcript_joiner(g_, long_length);
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

    void OutputPaths(const PathContainer& paths, const string& filename_base) const {
        WritePathsToFASTA(paths, filename_base);
    }


};


class PathInfoWriter {
protected:
    DECL_LOGGER("PathExtendIO")

public:

    void WritePaths(const PathContainer &paths, const string &filename) const {
        std::ofstream oss(filename.c_str());

        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            iter.get()->Print(oss);
        }

        oss.close();
    }
};

}

#endif /* PE_IO_HPP_ */
