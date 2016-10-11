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
    double coverage_;
    IOContigStorage(const std::string& sequence, double coverage) :
    sequence_(sequence), coverage_(coverage) { }
};

struct IOContigStorageGreater
{
    bool operator()(IOContigStorage const &a, IOContigStorage const &b) const {
        return a.sequence_.length() > b.sequence_.length();
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
    bool plasmid_contig_naming_;

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
                        for (size_t j = 0; j < abs(overlapLen); ++j) {
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
                 bool plasmid_contig_naming = false):
        g_(g), constructor_(constructor), k_(g.k()),
        ids_(), c_counter_(c_counter),
        plasmid_contig_naming_(plasmid_contig_naming)
    {
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
                           bool write_fastg = true) const {

        INFO("Writing contigs to " << filename_base);
        io::osequencestream_simple oss(filename_base + ".fasta");

        std::ofstream os_fastg;
        if (write_fastg)
            os_fastg.open((filename_base + ".paths").c_str());
        vector<IOContigStorage> storage;
        int i = 0;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (iter.get()->Length() <= 0)
                continue;
            i++;
            DEBUG("NODE " << i);
            BidirectionalPath* path = iter.get();
            path->Print();
            string contig_id;
            string path_string = ToString(*path);
            if (plasmid_contig_naming_) {
                EdgeId e = path->At(0);
                size_t component = c_counter_.GetComponent(e);
                contig_id = io::MakeContigComponentId(i, path_string.length(), path->Coverage(), component);
                oss.set_header(contig_id);
                oss << path_string;
            } else {
                if (path_string.length() >= g_.k()) {
                    IOContigStorage precontig(path_string, path->Coverage());
                    storage.push_back(precontig);
                }
            }

            if (write_fastg) {
                os_fastg << contig_id<< endl;
                os_fastg << ToFASTGString(*iter.get()) << endl;
                os_fastg << contig_id << "'" << endl;
                os_fastg << ToFASTGString(*iter.getConjugate()) << endl;
            }

        }
        std::sort(storage.begin(), storage.end(), IOContigStorageGreater());
        i = 0;
        for (auto precontig : storage) {
            ++i;
            std::string contig_id = io::MakeContigId(i, precontig.sequence_.length(), precontig.coverage_);
            oss.set_header(contig_id);
            oss << precontig.sequence_;
        }
        if (write_fastg)
            os_fastg.close();
        DEBUG("Contigs written");
    }

    void WriteFASTGPaths(const PathContainer& paths, const string& filename) const {
        INFO("Writing FASTG paths to " << filename);
        std::ofstream oss(filename.c_str());
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (iter.get()->Length() <= 0)
                continue;
            oss << ToFASTGString(*iter.get()) << endl;
            oss << ToFASTGString(*iter.getConjugate()) << endl;
        }
        oss.close();
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
