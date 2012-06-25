#pragma once
#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"
#include "standard_vis.hpp"

namespace online_visualization {

    struct Environment {
        private :
            const string name;
            size_t picture_counter_;
            string folder_;
            string file_name_base_;
            size_t max_vertices_;
            size_t edge_length_bound_;

            GraphPack gp_;
            map<EdgeId, string> coloring_;
            EdgePos positions_;
            CompositeLabeler<Graph> labeler_;

        public :
            
            Environment(const string& env_name) :
                name(env_name), 
                picture_counter_(0), 
                folder_("pictures"), 
                file_name_base_("picture"), 
                max_vertices_(40), 
                edge_length_bound_(1000), 
            {
                
            }
            
            void LoadFromGP(const GraphPack& gp) {
                gp_ = gp;
                //Loading Genome
                NewPathColorer<Graph> colorer(gp.g);
                NewExtendedSequenceMapper<K + 1, Graph> mapper(gp_.g, gp_.index, gp_.kmer_mapper);
                MappingPath<EdgeId> path1 = mapper.MapSequence(gp.genome);
                MappingPath<EdgeId> path2 = mapper.MapSequence(!gp.genome);
                coloring_ = colorer.ColorPath(path1.simple_path(), path2.simple_path());
                PosFiller<Graph, NewExtendedSequenceMapper<K + 1, Graph> > filler(gp.g, mapper, positions_);
                filler.Process(gp.genome, "ref0");
                filler.Process(!gp.genome, "ref1");
                   
            }

    };
}
