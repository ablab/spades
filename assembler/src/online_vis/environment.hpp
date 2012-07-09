#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "debruijn_stats.hpp"

namespace online_visualization {

    typedef debruijn_graph::NewExtendedSequenceMapper<debruijn_graph::K + 1, Graph> MapperClass;

    typedef debruijn_graph::PosFiller<Graph, MapperClass> FillerClass;

    typedef map<EdgeId, string> ColoringClass;

    struct Environment {
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
                gp_(cfg::get().ds.reference_genome, cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling),
                mapper_(gp_.g, gp_.index, gp_.kmer_mapper),
                filler_(gp_.g, mapper_, gp_.edge_pos),
                graph_struct_(gp_.g, &gp_.int_ids, &gp.edge_pos), 
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

            void ResetPositions() {
                mapper_(gp_.g, gp_.index, gp_.kmer_mapper);
                filler_(gp_.g, mapper_, gp_.edge_pos);
                filler_.Process(gp_.genome, "ref0");
                filler_.Process(!gp_.genome, "ref1");
            }

            string str() {
                stringstream ss;
                ss << name + " " + path;
                return ss.str();   
            }

            void set_max_vertices(size_t max_vertices) {
                max_vertices_ = max_vertices;
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
