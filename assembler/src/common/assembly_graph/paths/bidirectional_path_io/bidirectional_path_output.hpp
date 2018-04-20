//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/utils/edge_namer.hpp"
#include "io/graph/gfa_writer.hpp"
#include "io_support.hpp"

namespace path_extend {

template<class Graph>
class FastgWriter {
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    io::CanonicalEdgeHelper<Graph> short_namer_;
    io::CanonicalEdgeHelper<Graph> extended_namer_;

    std::string ToPathString(const BidirectionalPath &path) const {
        if (path.Empty())
            return "";
        std::string res = short_namer_.EdgeOrientationString(path.Front());
        for (size_t i = 1; i < path.Size(); ++i) {
            if (graph_.EdgeEnd(path[i - 1]) != graph_.EdgeStart(path[i]) || path.GapAt(i).gap > 0) {
                res += ";\n" + short_namer_.EdgeOrientationString(path[i]);
            } else {
                res += "," + short_namer_.EdgeOrientationString(path[i]);
            }
        }
        return res;
    }

    std::string FormHeader(const std::string &id,
                           const std::set<std::string>& next_ids) {
        std::stringstream ss;
        ss << id;
        if (next_ids.size() > 0) {
            auto delim = ":";
            for (const auto &s : next_ids) {
                ss  << delim << s;
                delim = ",";
            }
        }
        ss << ";";
        return ss.str();
    }

public:
    FastgWriter(const Graph &graph,
                io::EdgeNamingF<Graph> edge_naming_f = io::BasicNamingF<Graph>())
            : graph_(graph),
              short_namer_(graph_),
              extended_namer_(graph_, edge_naming_f, "", "'") {
    }

    void WriteSegmentsAndLinks(const std::string &fn) {
        io::OFastaReadStream os(fn);
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            std::set<std::string> next;
            for (EdgeId next_e : graph_.OutgoingEdges(graph_.EdgeEnd(e))) {
                next.insert(extended_namer_.EdgeOrientationString(next_e));
            }
            os << io::SingleRead(FormHeader(extended_namer_.EdgeOrientationString(e), next),
                                 graph_.EdgeNucls(e).str());
        }
    }

    void WritePaths(const ScaffoldStorage &scaffold_storage, const std::string &fn) const {
        std::ofstream os(fn);
        for (const auto& scaffold_info : scaffold_storage) {
            os << scaffold_info.name << "\n";
            os << ToPathString(*scaffold_info.path) << "\n";
            os << scaffold_info.name << "'" << "\n";
            os << ToPathString(*scaffold_info.path->GetConjPath()) << "\n";
        }
    }

};

class GFAPathWriter : public gfa::GFAWriter {
    void WritePath(const std::string& name, size_t segment_id,
                   const std::vector<std::string> &edge_strs) {
        os_ << "P" << "\t" ;
        os_ << name << "_" << segment_id << "\t";
        std::string delimeter = "";
        for (const auto& e : edge_strs) {
            os_ << delimeter << e;
            delimeter = ",";
        }
        os_ << "\t*\n";
//        delimeter = "";
//        for (size_t i = 0; i < edge_strs.size() - 1; ++i) {
//            os_ << delimeter << "*";
//            delimeter = ",";
//        }
//        os_ << "\n";
    }

public:
    using gfa::GFAWriter::GFAWriter;

    void WritePaths(const ScaffoldStorage &scaffold_storage) {
        for (const auto& scaffold_info : scaffold_storage) {
            const path_extend::BidirectionalPath &p = *scaffold_info.path;
            if (p.Size() == 0) {
                continue;
            }
            std::vector<std::string> segmented_path;
            //size_t id = p.GetId();
            size_t segment_id = 1;
            for (size_t i = 0; i < p.Size() - 1; ++i) {
                EdgeId e = p[i];
                segmented_path.push_back(edge_namer_.EdgeOrientationString(e));
                if (graph_.EdgeEnd(e) != graph_.EdgeStart(p[i+1]) || p.GapAt(i+1).gap > 0) {
                    WritePath(scaffold_info.name, segment_id, segmented_path);
                    segment_id++;
                    segmented_path.clear();
                }
            }

            segmented_path.push_back(edge_namer_.EdgeOrientationString(p.Back()));
            WritePath(scaffold_info.name, segment_id, segmented_path);
        }
    }

};

typedef std::function<void (const ScaffoldStorage&)> PathsWriterT;

class ContigWriter {
    const Graph& g_;
    shared_ptr<ContigNameGenerator> name_generator_;

public:
    static void WriteScaffolds(const ScaffoldStorage &scaffold_storage, const std::string &fn) {
        io::OFastaReadStream oss(fn);
        std::ofstream os_fastg;

        for (const auto& scaffold_info : scaffold_storage) {
            TRACE("Scaffold " << scaffold_info.name << " originates from path " << scaffold_info.path->str());
            oss << io::SingleRead(scaffold_info.name, scaffold_info.sequence);
        }
    }

    static PathsWriterT BasicFastaWriter(const std::string &fn) {
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

    void OutputPaths(const PathContainer &paths, const std::string &fn) const {
        OutputPaths(paths, BasicFastaWriter(fn));
    }

private:
    DECL_LOGGER("ContigWriter")
};

}
