//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_support.hpp"

namespace path_extend {
template<class Graph>
class FastgWriter {
    typedef typename Graph::EdgeId EdgeId;

    struct ExtendedContigId {
        string full_id_;
        string short_id_;

        ExtendedContigId() {}

        ExtendedContigId(string full_id, string short_id):
                full_id_(full_id), short_id_(short_id) {}
    };

    const Graph &graph_;
    map<EdgeId, ExtendedContigId> ids_;

    void FillEdgeIdMap(const ConnectedComponentCounter &cc_counter,
                                                const string &prefix) {
        int counter = 0;
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            if (ids_.count(e) == 0) {
                string id;
                if (cfg::get().pd) {
                    size_t c_id = cc_counter.GetComponent(e);
                    id = io::MakeContigComponentId(++counter, graph_.length(e) + graph_.k(), graph_.coverage(e), c_id, prefix);
                } else {
                    id = io::MakeContigId(++counter, graph_.length(e) + graph_.k(), graph_.coverage(e), prefix);
                }
                ids_[e] = ExtendedContigId(id, std::to_string(counter) + "+");
                if (e != graph_.conjugate(e))
                    ids_[graph_.conjugate(e)] =  ExtendedContigId(id + "'", std::to_string(counter) + "-");
            }
        }
    }


    string ToPathString(const BidirectionalPath &path) const {
        if (path.Empty())
            return "";
        string res = ids_.at(path.Front()).short_id_;
        for (size_t i = 1; i < path.Size(); ++i) {
            if (graph_.EdgeEnd(path[i - 1]) != graph_.EdgeStart(path[i])) {
                res += ";\n" + ids_.at(path[i]).short_id_;
            } else {
                res += "," + ids_.at(path[i]).short_id_;
            }
        }
        return res;
    }

    void ReportEdge(io::osequencestream_for_fastg& oss,
                    const string& sequence,
                    const string& id,
                    const set<string>& nex_ids) {
        oss.set_header(id);
        oss << nex_ids;
        oss << sequence;
    }

public:
    FastgWriter(const Graph &graph,
              const ConnectedComponentCounter &cc_counter,
              const string &prefix = "NODE")
            : graph_(graph) {
        FillEdgeIdMap(cc_counter, prefix);
    }

    void WriteSegmentsAndLinks(const string &fn) {
        io::osequencestream_for_fastg os(fn);
        for (auto it = graph_.ConstEdgeBegin(true); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            set<string> next;
            for (EdgeId next_e : graph_.OutgoingEdges(graph_.EdgeEnd(e))) {
                next.insert(ids_[next_e].full_id_);
            }
            ReportEdge(os, graph_.EdgeNucls(e).str(), ids_[e].full_id_, next);
            if (e != graph_.conjugate(e)) {
                set<string> next_conj;
                for (EdgeId next_e : graph_.OutgoingEdges(graph_.EdgeEnd(graph_.conjugate(e)))) {
                    next_conj.insert(ids_[next_e].full_id_);
                }
                ReportEdge(os, graph_.EdgeNucls(graph_.conjugate(e)).str(), ids_[graph_.conjugate(e)].full_id_, next_conj);
            }
        }
    }

    void WritePaths(const ScaffoldStorage &scaffold_storage, const string &fn) const {
        std::ofstream os(fn);
        for (const auto& scaffold_info : scaffold_storage) {
            os << scaffold_info.name << endl;
            os << ToPathString(*scaffold_info.path) << endl;
            os << scaffold_info.name << "'" << endl;
            os << ToPathString(*scaffold_info.path->GetConjPath()) << endl;
        }
    }

};

void WriteScaffolds(const ScaffoldStorage &scaffold_storage, const string &fn);

typedef std::function<void (const ScaffoldStorage&)> PathsWriterT;

//INFO("Writing contigs to " << filename_base);
//
//WriteScaffolds(storage, filename_base + ".fasta");
//inline void OutputContigsToFASTG(const Graph& g,
//                                 const string& contigs_output_filename,
//                                 const ConnectedComponentCounter &cc_counter) {
//    INFO("Outputting graph to " << contigs_output_filename << ".fastg");
//    ContigPrinter<Graph>(g).PrintContigsFASTG(ossfg, cc_counter);
//}

class ContigWriter {
    const Graph& g_;
    shared_ptr<ContigNameGenerator> name_generator_;

public:

    static PathsWriterT BasicFastaWriter(const string &fn) {
        return [=](const ScaffoldStorage& scaffold_storage) {
            WriteScaffolds(scaffold_storage, fn);
        };
    }

    ContigWriter(const Graph& g,
                 shared_ptr<ContigNameGenerator> name_generator) :
            g_(g),
            name_generator_(name_generator) {
    }

    void OutputPaths(const PathContainer &paths, const vector<PathsWriterT>& writers) const;

    void OutputPaths(const PathContainer &paths, PathsWriterT writer) const {
        OutputPaths(paths, vector<PathsWriterT>{writer});
    }

    void OutputPaths(const PathContainer &paths, const string &fn) const {
        OutputPaths(paths, BasicFastaWriter(fn));
    }

private:
    DECL_LOGGER("PathExtendIO")
};

}
