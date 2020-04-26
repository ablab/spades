//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "environment.hpp"
#include "errors.hpp"
#include "vis_logger.hpp"
#include "io/binary/graph_pack.hpp"
#include "modules/alignment/kmer_mapper.hpp"
#include "sequence/genome_storage.hpp"
#include "pipeline/config_struct.hpp"

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
              element_finder_(gp_.get<Graph>()),
              mapper_(MapperInstance(gp_)),
              filler_(gp_.get<Graph>(), mapper_, gp_.get_mutable<EdgePos>()),
              labeler_(gp_.get<Graph>(), gp_.get<EdgePos>()),
              path_finder_(gp_.get<Graph>()) {
            DEBUG("Environment constructor");
            gp_.get_mutable<debruijn_graph::KmerMapper<Graph>>().Attach();
            io::binary::BasePackIO().Load(path_, gp_);
//            debruijn_graph::graphio::ScanGraphPack(path_, gp_);
            DEBUG("Graph pack created")
            LoadFromGP();
        }

        inline bool IsCorrect() const {
            if (!CheckFileExists(path_ + ".grseq"))
                return false;

            size_t K = gp_.k();
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
            auto &edge_pos = gp_.get_mutable<EdgesPositionHandler<Graph>>();

            if (!edge_pos.IsAttached()) {
                edge_pos.Attach();
            }

            //Loading Genome and Handlers
            DEBUG("Colorer done");
            const auto &genome = gp_.get<GenomeStorage>();
            Path<EdgeId> path1 = mapper_->MapSequence(genome.GetSequence()).path();
            Path<EdgeId> path2 = mapper_->MapSequence(!genome.GetSequence()).path();
            coloring_ = visualization::graph_colorer::DefaultColorer(gp_.get<Graph>(), path1, path2);
            ResetPositions();
        }

        void LoadNewGenome(const Sequence& genome) {
            gp_.get_mutable<GenomeStorage>().SetSequence(genome);
            ResetPositions();
        }

        void ResetPositions() {
            auto &edge_pos = gp_.get_mutable<EdgesPositionHandler<Graph>>();

            if (!edge_pos.IsAttached())
                edge_pos.Attach();

            edge_pos.clear();

            const auto &genome = gp_.get<GenomeStorage>();
            filler_.Process(genome.GetSequence(), "ref0");
            filler_.Process(!genome.GetSequence(), "ref1");
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
            return gp_.k();
        }

        const Graph& graph() const {
            return gp_.get<Graph>();
        }

        GraphPack& graph_pack() {
            return gp_;
        }

        Graph& graph() {
            return gp_.get_mutable<Graph>();;
        }

        Sequence genome() const {
            return gp_.get<GenomeStorage>().GetSequence();
        }

        const MapperClass& mapper() const {
            return *mapper_;
        }

        const debruijn_graph::ReadPathFinder<Graph>& path_finder() const {
                    return path_finder_;
        }

        const Index &index() const {
            return gp_.get<Index>();
        }

        const KmerMapperClass& kmer_mapper() const {
            return gp_.get<KmerMapperClass>();
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
