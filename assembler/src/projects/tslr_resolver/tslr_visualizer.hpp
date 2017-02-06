#pragma once

#include "common/barcode_index/barcode_info_extractor.hpp"

namespace barcode_index {

    template<class G>
    class BarcodeDistGraphLabeler : public visualization::graph_labeler::StrGraphLabeler<Graph> {
        typedef visualization::graph_labeler::StrGraphLabeler <Graph> base;
        typedef typename G::EdgeId EdgeId;
        typedef typename G::VertexId VertexId;
        typedef AbstractBarcodeIndex b_mapper;
        typedef std::set <EdgeId> edge_set_t;
    private:
        const shared_ptr <b_mapper> barcode_mapper_;
        //Code duplication to avoid passing unique edge storage to Labeler.

        bool ConservativeByTopology(EdgeId e, const Graph& g ) const {
            size_t incoming = g.IncomingEdgeCount(g.EdgeStart(e));
            size_t outcoming = g.OutgoingEdgeCount(g.EdgeEnd(e));
            return incoming < 2 and outcoming < 2;
        }

        bool ConservativeByLength(EdgeId e, const Graph& g, size_t threshold) const {
            return g.length(e) > threshold;
        }
    public:
        BarcodeDistGraphLabeler(const G &g, shared_ptr <b_mapper>& mapper) :
                base(g), barcode_mapper_(mapper){}

        size_t barcode_threshold = 5;

        std::string label(EdgeId e) const {
            std::string ans;
            ans += ("Id: " + std::to_string(this->graph().int_id(e)) + ' '
                    +  "Length: " + std::to_string(this->graph().length(e)) + ' '
                    + "Coverage: " + std::to_string(this->graph().coverage(e)) + '\n');

            ans += ("Barcodes at head: " + std::to_string(barcode_mapper_->GetHeadBarcodeNumber(e)) + '\n');
            ans += "Barcodes at tail: " + std::to_string(barcode_mapper_->GetTailBarcodeNumber(e)) + '\n';

            bool is_unique = ConservativeByTopology(e, this->graph()) and
                    ConservativeByLength(e, this->graph(), 500);
            if(is_unique) {
                ans += "Unique\n";
            }
            else
                ans += "Non-unique\n";
            return ans;
        }
    };

    template <class Graph>
    class TslrVisualizer {
        typedef AbstractBarcodeIndex b_mapper;
    private:
        const debruijn_graph::conj_graph_pack &gp_;
        shared_ptr<b_mapper> barcode_mapper_;
    public:
        TslrVisualizer(const debruijn_graph::conj_graph_pack& gp, shared_ptr<b_mapper> barcode_mapper) :
                gp_(gp), barcode_mapper_(barcode_mapper) { }

        std::string pics_folder = cfg::get().output_dir + '/' + "pictures";

        void DrawEdgeComponent(const EdgeId &edge, size_t num) {
            size_t max_vertices = 60;
            size_t edge_length_bound = 10000;
            auto component = omnigraph::EdgeNeighborhood(gp_.g, edge, max_vertices, edge_length_bound);
            auto edge_set = component.edges();
            visualization::graph_labeler::EdgePosGraphLabeler <Graph> pos_labeler(gp_.g, gp_.edge_pos);
            BarcodeDistGraphLabeler<Graph> barcode_labeler(gp_.g, barcode_mapper_);
            visualization::graph_labeler::CompositeLabeler <Graph> composite_labeler(pos_labeler, barcode_labeler);
            auto colorer = visualization::graph_colorer::DefaultColorer(gp_.g);
            visualization::visualization_utils::WriteComponent(component,
                                                     pics_folder + "/edge_component_" + std::to_string(num) + "_.dot",
                                                     colorer, composite_labeler);
        }


        void DrawRandomRepeats() {
            srand(3);
            std::string test_file = pics_folder + "/random_repeats";
            std::ofstream fout;
            fout.open(test_file);
            size_t components_number = 20;
            std::vector <EdgeId> repeats;
            edge_it_helper edge_it(gp_.g);
            for (auto it = edge_it.begin(); it != edge_it.end(); ++it) {
                if (gp_.edge_pos.GetEdgePositions(*it).size() >= 2) {
                    repeats.push_back(*it);
                }
            }
            fout << repeats.size() << std::endl;
            for (auto e : repeats) {
                fout << gp_.g.str(e) << "  " << gp_.edge_pos.GetEdgePositions(e).size() << std::endl;
            }
            size_t size = repeats.size();
            VERIFY(size > 1);
            for (size_t i = 0; i < components_number; ++i) {
                size_t ind = rand() % size;
                auto it = repeats.begin();
                for (size_t j = 0; j < ind; ++j) {
                    it++;
                }
                DrawEdgeComponent(*it, i);
            }
        }

    };
}
