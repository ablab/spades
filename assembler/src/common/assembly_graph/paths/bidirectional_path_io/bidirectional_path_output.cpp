//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bidirectional_path_output.hpp"
#include "assembly_graph/core/graph.hpp"
#include "io/graph/gfa_writer.hpp"

#include <unordered_set>

namespace path_extend {

std::string PathWriter::ToPathString(const BidirectionalPath &path) const {
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

void FastgPathWriter::WritePaths(const ScaffoldStorage &scaffold_storage, const std::filesystem::path &fn) const {
    std::ofstream os(fn);
    for (const auto& scaffold_info : scaffold_storage) {
        os << scaffold_info.name << "\n"
           << path_writer_.ToPathString(*scaffold_info.path) << "\n"
           << scaffold_info.name << "'" << "\n"
           << path_writer_.ToPathString(*scaffold_info.path->GetConjPath()) << "\n";
    }
}

void GFAPathWriter::WritePath(const std::string &name, size_t segment_id,
                              const std::vector<std::string> &edge_strs,
                              const std::string &flags) {
    os_ << 'P' << '\t' ;
    os_ << name << '_' << segment_id << '\t';
    std::string delimeter = "";
    for (const auto& e : edge_strs) {
        os_ << delimeter << e;
        delimeter = ',';
    }
    os_ << "\t*";
    if (flags.length())
        os_ << '\t' << flags;
    os_ << '\n';
}

void GFAPathWriter::WritePath(const std::string &name,
                              const std::vector<std::string> &edge_strs,
                              const std::string &flags) {
    os_ << 'P' << '\t' ;
    os_ << name << '\t';
    std::string delimeter = "";
    for (const auto& e : edge_strs) {
        os_ << delimeter << e;
        delimeter = ',';
    }
    os_ << "\t*";
    if (flags.length())
        os_ << '\t' << flags;
    os_ << '\n';
}

GFAPathWriter::GFAPathWriter(const Graph &graph, std::ostream &os,
                             io::EdgeNamingF<Graph> naming_f,
                             Version version)
        : gfa::GFAWriter(graph, os, naming_f),
          version_(version) {}

void GFAPathWriter::WritePaths11(const std::vector<EdgeId> &edges,
                                 const std::string &name,
                                 const std::string &flags) {
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

void GFAPathWriter::WriteJumpLinks(const JumpLinks &jump_links) {
    std::vector<std::pair<EdgeId, EdgeId>> sorted_jump_links;
    std::copy(jump_links.begin(), jump_links.end(), std::back_inserter(sorted_jump_links));
    std::sort(sorted_jump_links.begin(), sorted_jump_links.end());
    for (auto &link : sorted_jump_links) {
        os_ << 'J' << '\t'
            << edge_namer_.EdgeOrientationString(link.first, "\t") << '\t'
            << edge_namer_.EdgeOrientationString(link.second, "\t") << '\t'
            << "*" << '\t'
            << "SC:i:1"
            << '\n';
    }
}


void GFAPathWriter::WritePaths12(const std::vector<EdgeId> &edges,
                                 const std::string &name,
                                 const std::string &flags) {
    JumpLinks jump_links;

    // First, determine the set of jump links to emit
    for (size_t i = 0; i < edges.size() - 1; ++i) {
        EdgeId e = edges[i];
        if (graph_.EdgeEnd(e) != graph_.EdgeStart(edges[i+1]))
            jump_links.emplace(edges[i], edges[i+1]);
    }

    // Emit jump links
    WriteJumpLinks(jump_links);

    // Emit path. We cannot use function above as we're single-segmented
    os_ << 'P' << '\t'
        << name << '\t';
    std::string delimiter = "";
    for (size_t i = 0; i < edges.size() - 1; ++i) {
        EdgeId e = edges[i];
        os_ << delimiter << edge_namer_.EdgeOrientationString(e);
        delimiter = (graph_.EdgeEnd(e) == graph_.EdgeStart(edges[i+1]) ? "," : ";");
    }
    os_ << delimiter << edge_namer_.EdgeOrientationString(edges.back());

    os_ << "\t*";
    if (flags.length())
        os_ << '\t' << flags;
    os_ << '\n';
}

void GFAPathWriter::WritePaths(const std::vector<EdgeId> &edges,
                               const std::string &name,
                               const std::string &flags) {
    if (version_ == Version::GFAv11)
        WritePaths11(edges, name, flags);
    else
        WritePaths12(edges, name, flags);;
}


void GFAPathWriter::WritePaths(const ScaffoldStorage &scaffold_storage) {
    if (version_ == Version::GFAv11)
        WritePaths11(scaffold_storage);
    else
        WritePaths12(scaffold_storage);
}

void GFAPathWriter::WritePaths(const gfa::GFAReader::GFAPath &path,
                               const std::string &flags) {
    std::vector<std::string> path_segment;

    for (EdgeId e : path.edges)
        path_segment.push_back(this->edge_namer_.EdgeOrientationString(e));

    WritePath(path.name, path_segment, flags);
}

void GFAPathWriter::WritePaths12(const ScaffoldStorage &scaffold_storage) {
    JumpLinks jump_links;

    for (const auto& scaffold_info : scaffold_storage) {
        const path_extend::BidirectionalPath &p = *scaffold_info.path;
        if (p.Size() == 0)
            continue;

        for (size_t i = 0; i < p.Size() - 1; ++i) {
            EdgeId e = p[i];
            if (graph_.EdgeEnd(e) != graph_.EdgeStart(p[i+1]) || p.GapAt(i+1).gap > 0)
                jump_links.emplace(e, p[i+1]);
        }
    }

    WriteJumpLinks(jump_links);

    for (const auto& scaffold_info : scaffold_storage) {
        const path_extend::BidirectionalPath &p = *scaffold_info.path;
        if (p.Size() == 0)
            continue;

        os_ << 'P' << '\t'
            << scaffold_info.name << '\t';
        std::string delimiter = "";
        for (size_t i = 0; i < p.Size() - 1; ++i) {
            EdgeId e = p[i];
            os_ << delimiter << edge_namer_.EdgeOrientationString(e);
            delimiter = (graph_.EdgeEnd(e) == graph_.EdgeStart(p[i+1]) ? "," : ";");
        }
        os_ << delimiter << edge_namer_.EdgeOrientationString(p.Back())
            << '\t' << '*';
        if (p.IsCircular())
            os_ << '\t' << "TP:Z:circular";
        os_ << '\n';
    }
}

void GFAPathWriter::WritePaths11(const ScaffoldStorage &scaffold_storage) {
    for (const auto& scaffold_info : scaffold_storage) {
        const path_extend::BidirectionalPath &p = *scaffold_info.path;
        if (p.Size() == 0)
            continue;

        std::vector<std::string> segmented_path;
        // size_t id = p.GetId();
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
        WritePath(scaffold_info.name, segment_id, segmented_path,
                  segment_id == 1 && p.IsCircular() ? "TP:Z:circular" : "");
    }
}


void ContigWriter::OutputPaths(const PathContainer &paths, const std::vector<PathsWriterT> &writers) const {
    ScaffoldStorage storage;

    ScaffoldSequenceMaker scaffold_maker(g_);
    DEBUG("started" << paths.size());
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        const BidirectionalPath &path = iter.get();
        DEBUG("path: " <<  path.Length());
        if (path.Length() <= 0)
            continue;
        auto path_string = scaffold_maker.MakeSequence(path);
        if (path_string.length() >= g_.k()) {
            storage.emplace_back(path_string, &path);
        }
        DEBUG("over");
    }
    DEBUG("sort");
    //sorting by length and coverage
    std::sort(storage.begin(), storage.end(), [] (const ScaffoldInfo &a, const ScaffoldInfo &b) {
        if (a.length() == b.length())
            return math::gr(a.coverage(), b.coverage());
        return a.length() > b.length();
    });
    DEBUG("preprocess");
    name_generator_->Preprocess(paths);
    DEBUG(storage.size());
    for (size_t i = 0; i < storage.size(); ++i) {
        DEBUG("name");
        storage[i].name = name_generator_->MakeContigName(i+1, storage[i]);
    }
    DEBUG("wrt");
    for (auto& writer : writers) {
        writer(storage);
    }
    DEBUG("Contigs written");
}

}
