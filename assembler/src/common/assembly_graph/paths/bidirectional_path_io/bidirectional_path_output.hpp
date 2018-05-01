//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/utils/edge_namer.hpp"
#include "io/graph/gfa_writer.hpp"
#include "io/graph/fastg_writer.hpp"
#include "io_support.hpp"

namespace path_extend {

class FastgPathWriter : public io::FastgWriter {
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

public:
    using io::FastgWriter::FastgWriter;

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
                   const std::vector<std::string> &edge_strs,
                   const std::string &flags) {
        os_ << "P" << "\t" ;
        os_ << name << "_" << segment_id << "\t";
        std::string delimeter = "";
        for (const auto& e : edge_strs) {
            os_ << delimeter << e;
            delimeter = ",";
        }
        os_ << "\t*";
        if (flags.length())
            os_ << "\t" << flags;
        os_ << "\n";
    }

public:
    using gfa::GFAWriter::GFAWriter;

    void WritePaths(const std::vector<EdgeId> &edges,
                    const std::string &name, const std::string &flags = "") {
        std::vector<std::string> segmented_path;
        size_t segment_id = 1;
        for (size_t i = 0; i < edges.size() - 1; ++i) {
            EdgeId e = edges[i];
            segmented_path.push_back(edge_namer_.EdgeOrientationString(e));
            if (graph_.EdgeEnd(e) != graph_.EdgeStart(edges[i+1])) {
                WritePath(name, segment_id, segmented_path, flags);
                segment_id++;
                segmented_path.clear();
            }
        }

        segmented_path.push_back(edge_namer_.EdgeOrientationString(edges.back()));
        WritePath(name, segment_id, segmented_path, flags);
    }

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
                    WritePath(scaffold_info.name, segment_id, segmented_path, "");
                    segment_id++;
                    segmented_path.clear();
                }
            }

            segmented_path.push_back(edge_namer_.EdgeOrientationString(p.Back()));
            WritePath(scaffold_info.name, segment_id, segmented_path, "");
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
