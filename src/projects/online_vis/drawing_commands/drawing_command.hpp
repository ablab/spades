//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "projects/online_vis/argument_list.hpp"
#include "projects/online_vis/command.hpp"
#include "projects/online_vis/environment.hpp"
#include "projects/online_vis/errors.hpp"

#include "alignment/sequence_mapper.hpp"
#include "io/reads/single_read.hpp"

namespace online_visualization {

class DrawingCommand : public LocalCommand<DebruijnEnvironment> {
protected:
    void DrawPicture(DebruijnEnvironment& curr_env, VertexId vertex, string label = "") const {
        create_directory(curr_env.folder_);
        filesystem::path file_name = curr_env.folder_ / (curr_env.GetFormattedPictureCounter() + "_" + curr_env.file_name_base_ + "_" + label + "_" + ".dot");
        omnigraph::GraphComponent<Graph> component = VertexNeighborhood(curr_env.graph(), vertex, curr_env.max_vertices_, curr_env.edge_length_bound_);
        visualization::visualization_utils::WriteComponent<Graph>(component, file_name, curr_env.coloring_, curr_env.labeler());
        LOG("The picture is written to " << file_name);

        curr_env.picture_counter_++;
    }

    void DrawPicturesAlongPath(DebruijnEnvironment& curr_env, const vector<EdgeId>& path, string label = "", size_t length = 0) const {
        create_directory(curr_env.folder_);
        filesystem::path directory = curr_env.folder_ / (curr_env.GetFormattedPictureCounter() + "_" + "edges_" + to_string(path.size()) + "_" + "length_" + to_string(length) + "_" + curr_env.file_name_base_);
        create_directory(directory);
        filesystem::path directory_path = directory.concat(label + "_");
        visualization::visualization_utils::WriteComponentsAlongPath<Graph>(curr_env.graph(), path, directory_path, curr_env.coloring_, curr_env.labeler());
        LOG("The pictures is written to " << directory);

        curr_env.picture_counter_++;
    }

    void DrawPicturesAlongContig(DebruijnEnvironment& curr_env, io::SingleRead contig) const {
        string label = contig.name();
        DrawPicturesAlongPath(curr_env, curr_env.mapper().MapRead(contig).simple_path(), label, contig.size());
        LOG("Contig " << contig.name() << " has been drawn");
    }

    void DrawPicturesBetweenTwoEdges(DebruijnEnvironment& curr_env,
                                     const EdgeId& edge1, const EdgeId& edge2,
                                     const vector<EdgeId>& path, const size_t barcode_threshold,
                                     const string& label = "") const {
        size_t pos1 = GetIndexInPath(path, edge1);
        size_t pos2 = GetIndexInPath(path, edge2);
        if (pos1 == path.size() or pos2 == path.size()) {
            LOG("One of the given edges was not found in the contig.");
            return;
        }

        make_dir(curr_env.folder_);
        stringstream namestream;
        namestream << curr_env.folder_ << "/" << label << "_" << curr_env.file_name_base_  << "_" << curr_env.GetFormattedPictureCounter() << "/";
        string directory = namestream.str();
        make_dir(directory);
        namestream << label << "_";

        barcode_index::FrameBarcodeIndexInfoExtractor extractor(curr_env.gp_.barcode_mapper_ptr, curr_env.graph());
        vector<barcode_index::BarcodeId> shared_barcodes = extractor.GetSharedBarcodes(edge1, edge2);
        auto barcode_colorer_ptr =
                make_shared<visualization::graph_colorer::BarcodeCollectionColorer<Graph>>(shared_barcodes,
                                                                                    extractor,
                                                                                    barcode_threshold, "blue");
        auto colorer = curr_env.coloring_;
        auto edge_colorer = make_shared<visualization::graph_colorer::CompositeEdgeColorer<Graph>>("black");
        edge_colorer->AddColorer(colorer);
        edge_colorer->AddColorer(barcode_colorer_ptr);
        vector<EdgeId> arguments({edge1, edge2});
//        edge_colorer->AddColorer(make_shared<visualization::graph_colorer::SetColorer<Graph>>(curr_env.gp_.g,
//                                                                                              arguments, "darkgreen"));

        shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> barcode_colorer =
                make_shared<visualization::graph_colorer::CompositeGraphColorer<Graph>>(colorer, edge_colorer);

        visualization::visualization_utils::WriteComponentsAlongPath<Graph>(curr_env.graph(), GetSubPath(path, pos1, pos2),
                                                                            namestream.str(), barcode_colorer, curr_env.labeler());
        LOG("Gap within contig " << label << " between edges "
                                 << edge1.int_id() << " and " << edge2.int_id() << " has been drawn.");
        LOG("The pictures is written to " << directory);

        curr_env.picture_counter_++;
    }

    size_t GetIndexInPath(const vector<EdgeId>& path, const EdgeId& edge) const {
        for (size_t i = 0; i < path.size(); ++i) {
            if (path[i] == edge) return i;
        }
        return path.size();
    }

    const vector<EdgeId> GetSubPath(const vector<EdgeId>& path, size_t pos1, size_t pos2) const {
        vector<EdgeId> result;
        VERIFY(pos1 < path.size() and pos2 < path.size());
        for (size_t i = std::min(pos1, pos2); i <= std::max(pos1, pos2); ++i) {
            result.push_back(path[i]);
        }
        return result;
    }

    void DrawPicturesInContigGap(DebruijnEnvironment& curr_env, io::SingleRead contig,
                           const EdgeId& edge1, const EdgeId& edge2, const size_t barcode_threshold) const {
        string label = contig.name();
        DrawPicturesBetweenTwoEdges(curr_env, edge1, edge2,
                                    curr_env.mapper().MapRead(contig).simple_path(), barcode_threshold, label);
    }

    void DrawConnectedComponents (DebruijnEnvironment& curr_env,  int min_size, int max_size, string label = "") const {
        create_directory(curr_env.folder_);
        filesystem::path name = curr_env.folder_ / (curr_env.GetFormattedPictureCounter() + "_" + curr_env.file_name_base_);
        create_directory(name);
        name.concat(label);
        create_directory(name);
        visualization::visualization_utils::WriteSizeLimitedComponents<Graph>(curr_env.graph(), name, omnigraph::ConnectedSplitter<Graph>(curr_env.graph()), curr_env.coloring_, curr_env.labeler(), min_size, max_size, 10000000);
        LOG("The pictures is written to " << name);
        curr_env.picture_counter_++;
    }

    //TODO: copy zgrviewer
    int ShowPicture(DebruijnEnvironment& curr_env, VertexId vertex, string label = "") const {
        DrawPicture(curr_env, vertex, label);
        stringstream command_line_string;
        command_line_string << "gnome-open " << curr_env.folder_ << "/" << curr_env.file_name_base_
                            << "_" << label << "_" << curr_env.GetFormattedPictureCounter()
                            << "_*_.dot & > /dev/null < /dev/null";
        int result = system(command_line_string.str().c_str());

        return result;
    }

    void DrawVertex(DebruijnEnvironment& curr_env, size_t vertex_id, string label = "") const {
        DrawPicture(curr_env, curr_env.finder().ReturnVertexId(vertex_id), label);
    }


public:
    DrawingCommand(string command_type) : LocalCommand<DebruijnEnvironment>(command_type)
    {
    }

    virtual ~DrawingCommand()
    {
    }
};
}
