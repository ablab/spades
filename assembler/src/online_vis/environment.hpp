#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"
#include "standard.hpp"

namespace online_visualization {

    struct Environment {
        typedef debruijn_graph::NewExtendedSequenceMapper<debruijn_graph::K + 1, Graph> MapperClass;
        typedef debruijn_graph::PosFiller<Graph, MapperClass> FillerClass;
        private :
            const string name;
            size_t picture_counter_;
            string folder_;
            string file_name_base_;
            size_t max_vertices_;
            size_t edge_length_bound_;

            map<EdgeId, string> coloring_;
            GraphPack gp_;

            //EdgePos positions_;
            //CompositeLabeler<Graph> labeler_;

        public :
            Environment(const string& env_name, const string& path) :
                name(env_name), 
                picture_counter_(0), 
                folder_("pictures_" + env_name), 
                file_name_base_("picture_" + env_name), 
                max_vertices_(40), 
                edge_length_bound_(1000)
            {
                Sequence genome = cfg::get().ds.reference_genome;
                GraphPack gp_(genome, cfg::get().pos.max_single_gap, cfg::get().pos.careful_labeling);
                ScanGraphPack(path, gp_);
                LoadFromGP();
            }
            
            void LoadFromGP() {
                //Loading Genome and Handlers
                NewPathColorer<Graph> colorer(gp_.g);
                MapperClass mapper(gp_.g, gp_.index, gp_.kmer_mapper);
                MappingPath<EdgeId> path1 = mapper.MapSequence(gp_.genome);
                MappingPath<EdgeId> path2 = mapper.MapSequence(!gp_.genome);
                coloring_ = colorer.ColorPath(path1.simple_path(), path2.simple_path());
                FillerClass filler(gp_.g, mapper, gp_.edge_pos);
                filler.Process(gp_.genome, "ref0");
                filler.Process(!gp_.genome, "ref1");
                   
            }

    };
}
