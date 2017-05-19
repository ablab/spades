//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "environment.hpp"
#include "pipeline/graphio.hpp"
namespace online_visualization {

class DebruijnEnvironment : public Environment {
    friend class DrawingCommand;

    private :
        size_t picture_counter_;
        string folder_;
        string file_name_base_;
        size_t max_vertices_;
        size_t edge_length_bound_;

        GraphPack gp_;
        GraphElementFinder<Graph> element_finder_;
        std::shared_ptr<MapperClass> mapper_;
        FillerClass filler_;
        visualization::graph_labeler::DefaultLabeler<Graph> labeler_;
        debruijn_graph::ReadPathFinder<Graph> path_finder_;
        ColoringClass coloring_;
        //CompositeLabeler<Graph> labeler_;

    public :

        typedef debruijn_graph::Index EdgeIndexT;

        DebruijnEnvironment(const string& env_name, const string& env_path, size_t K = cfg::get().K)
            : Environment(env_name, env_path),
              picture_counter_(0),
              folder_("pictures_" + name_),
              file_name_base_("picture"),
              max_vertices_(40),
              edge_length_bound_(1000),
              gp_(K, "./tmp", cfg::get().ds.reads.lib_count(), 
                  std::vector<std::string>(0),
                  cfg::get().flanking_range,
                  cfg::get().pos.max_mapping_gap,
                  cfg::get().pos.max_gap_diff),
              element_finder_(gp_.g),
              mapper_(new MapperClass(gp_.g, gp_.index, gp_.kmer_mapper)),
              filler_(gp_.g, mapper_, gp_.edge_pos),
              labeler_(gp_.g, gp_.edge_pos),
              path_finder_(gp_.g) {
            DEBUG("Environment constructor");
            gp_.kmer_mapper.Attach();
            debruijn_graph::graphio::ScanGraphPack(path_, gp_);
//            debruijn_graph::graphio::ScanGraphPack(path_, gp_);
            DEBUG("Graph pack created")
            LoadFromGP();
        }

        inline bool IsCorrect() const {
            if (!CheckFileExists(path_ + ".grp"))
                return false;
            if (!CheckFileExists(path_ + ".sqn"))
                return false;

            size_t K = gp_.k_value;
            if (!(K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K)) {
                LOG("K " << K << " is out of bounds");
                    return false;
            }
            if (K % 2 == 0) {
                LOG("K must be odd");
                    return false;
            }

            return true;
        }

        void LoadFromGP() {
            if (!gp_.edge_pos.IsAttached()) {
                gp_.edge_pos.Attach();
            }

            //Loading Genome and Handlers
            DEBUG("Colorer done");
            Path<EdgeId> path1 = mapper_->MapSequence(gp_.genome.GetSequence()).path();
            Path<EdgeId> path2 = mapper_->MapSequence(!gp_.genome.GetSequence()).path();
            coloring_ = visualization::graph_colorer::DefaultColorer(gp_.g, path1, path2);
            ResetPositions();
        }

        void LoadNewGenome(const Sequence& genome) {
            gp_.genome.SetSequence(genome);
            ResetPositions();
        }

        void ResetPositions() {
            if (!gp_.edge_pos.IsAttached())
                gp_.edge_pos.Attach();

            gp_.edge_pos.clear();
            filler_.Process(gp_.genome.GetSequence(), "ref0");
            filler_.Process(!gp_.genome.GetSequence(), "ref1");
        }

        string GetFormattedPictureCounter() const {
            stringstream tmpstream;
            size_t number_of_digs = 0;
            size_t pc = picture_counter_;

            do {
                pc /= 10;
                number_of_digs++;
            } while (pc > 0);

            for (size_t i = 0; i < 4 - number_of_digs; ++i)
                tmpstream << '0';
            tmpstream << picture_counter_;
            return tmpstream.str();
        }

        void inc_pic_counter() {
            picture_counter_++;
        }

        size_t k_value() const {
            return gp_.k_value;
        }

        const Graph& graph() const {
            return gp_.g;
        }

        GraphPack& graph_pack() {
            return gp_;
        }

        Graph& graph() {
            return gp_.g;
        }

        Sequence genome() const {
            return gp_.genome.GetSequence();
        }

        const MapperClass& mapper() const {
            return *mapper_;
        }

        const debruijn_graph::ReadPathFinder<Graph>& path_finder() const {
                    return path_finder_;
        }

        const EdgeIndexT& index() const {
            return gp_.index;
        }

        const KmerMapperClass& kmer_mapper() const {
            return gp_.kmer_mapper;
        }

        const ElementFinder& finder() const {
            return element_finder_;
        }

        void set_max_vertices(size_t max_vertices) {
            max_vertices_ = max_vertices;
        }

        void set_folder(string folder) {
            folder_ = folder;
        }

        string folder() const {
            return folder_;
        }

        void set_file_name(string file_name) {
            file_name_base_ = file_name;
        }

        string file_name() {
            return file_name_base_;
        }

        size_t edge_length_bound() const {
            return edge_length_bound_;
        }

        FillerClass& filler() {
            return filler_;
        }

        visualization::graph_labeler::GraphLabeler<Graph>& labeler() {
            return labeler_;
        }

        ColoringClass& coloring() {
            return coloring_;
        }

};

}
