#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "debruijn_stats.hpp"

namespace online_visualization {


    typedef debruijn_graph::NewExtendedSequenceMapper<Graph> MapperClass;

    typedef debruijn_graph::PosFiller<Graph, MapperClass> FillerClass;

    typedef debruijn_graph::KmerMapper<Graph> KmerMapperClass;

    typedef map<EdgeId, string> ColoringClass;


    struct Environment {
        friend class DrawingCommand;
        private :
            const string name;
            const string path;
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
            Environment(const string& env_name, const string& env_path) :
                name(env_name),
                path(env_path),
                picture_counter_(0),
                folder_("pictures_" + name),
                file_name_base_("picture"),
                max_vertices_(40),
                edge_length_bound_(1000),
                gp_(cfg::get().K, env_path + "/tmp", cfg::get().ds.reference_genome, cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling),
                mapper_(gp_.g, gp_.index, gp_.kmer_mapper, cfg::get().K + 1),
                filler_(gp_.g, mapper_, gp_.edge_pos),
                graph_struct_(gp_.g, &gp_.int_ids, &gp_.edge_pos), 
                tot_lab_(&graph_struct_)
            {
                ScanGraphPack(path, gp_);
                LoadFromGP();
            }
            
            void LoadFromGP() {
                //Loading Genome and Handlers
                NewPathColorer<Graph> colorer(gp_.g);
                MappingPath<EdgeId> path1 = mapper_.MapSequence(gp_.genome);
                MappingPath<EdgeId> path2 = mapper_.MapSequence(!gp_.genome);
                coloring_ = colorer.ColorPath(path1.simple_path(), path2.simple_path());
                ResetPositions();
            }

            void LoadNewGenome(const Sequence& genome) {
                gp_.genome = genome;
                ResetPositions();
            }

            void ResetPositions() {
                gp_.edge_pos.clear();
                MapperClass mapper_(gp_.g, gp_.index, gp_.kmer_mapper, gp_.k_value + 1);
                FillerClass filler_(gp_.g, mapper_, gp_.edge_pos);
                filler_.Process(gp_.genome, "ref0");
                filler_.Process(!gp_.genome, "ref1");
            }

            string str() const {
                stringstream ss;
                ss << name + " " + path;
                return ss.str();   
            }

            const Graph& graph() const {
                return gp_.g;
            }

            const Sequence& genome() const {
                return gp_.genome;   
            }

            const MapperClass& mapper() const {
                return mapper_;   
            }

            const debruijn_graph::EdgeIndex<Graph>& index() const {
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
