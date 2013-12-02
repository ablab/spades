//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "environment.hpp"

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
        MapperClass mapper_;
        FillerClass filler_;
        total_labeler_graph_struct graph_struct_;
        total_labeler tot_lab_;
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
              gp_(K, "./tmp", cfg::get().ds.reads.lib_count(), cfg::get().ds.reference_genome),
              mapper_(gp_.g, gp_.index, gp_.kmer_mapper),
              filler_(gp_.g, mapper_, gp_.edge_pos),
              graph_struct_(gp_.g, &gp_.int_ids, &gp_.edge_pos),
              tot_lab_(&graph_struct_) {

            DEBUG("Environment constructor");
            debruijn_graph::graphio::ScanGraphPack(path_, gp_);
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
            //Loading Genome and Handlers
            DEBUG("Colorer done");
            Path<EdgeId> path1 = mapper_.MapSequence(gp_.genome).path();
            Path<EdgeId> path2 = mapper_.MapSequence(!gp_.genome).path();
        	coloring_ = omnigraph::visualization::DefaultColorer(gp_.g, path1, path2);
            ResetPositions();
        }

        void LoadNewGenome(const Sequence& genome) {
            gp_.genome = genome;
            ResetPositions();
        }

        void ResetPositions() {
            gp_.edge_pos.clear();
            MapperClass mapper_(gp_.g, gp_.index, gp_.kmer_mapper);
            FillerClass filler_(gp_.g, mapper_, gp_.edge_pos);
            filler_.Process(gp_.genome, "ref0");
            filler_.Process(!gp_.genome, "ref1");
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

        size_t k_value() const {
            return gp_.k_value;
        }

        const Graph& graph() const {
            return gp_.g;
        }

        Graph& graph() {
            return gp_.g;
        }

        const Sequence& genome() const {
            return gp_.genome;
        }

        const MapperClass& mapper() const {
            return mapper_;
        }

        const EdgeIndexT& index() const {
            return gp_.index;
        }

        const KmerMapperClass& kmer_mapper() const {
            return gp_.kmer_mapper;
        }

        const IdTrackHandler<Graph>& int_ids() const {
            return gp_.int_ids;
        }

        void set_max_vertices(size_t max_vertices) {
            max_vertices_ = max_vertices;
        }

        void set_folder(string folder) {
            folder_ = folder;
        }

        void set_file_name(string file_name) {
            file_name_base_ = file_name;
        }

        FillerClass& filler() {
            return filler_;
        }

        total_labeler& tot_lab() {
            return tot_lab_;
        }

        ColoringClass& coloring() {
            return coloring_;
        }

};

}
