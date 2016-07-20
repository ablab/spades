#pragma once

#include "barcode_mapper.hpp"
#include "../../modules/visualization/graph_labeler.hpp"

namespace tslr_resolver {

    template<class G>
    class BarcodeDistGraphLabeler : public omnigraph::StrGraphLabeler<Graph> {
        typedef omnigraph::StrGraphLabeler <Graph> base;
        typedef typename G::EdgeId EdgeId;
        typedef typename G::VertexId VertexId;
        typedef BarcodeMapper b_mapper;
        typedef std::set <EdgeId> edge_set_t;
    private:
        const b_mapper &barcode_mapper_;
        const edge_set_t &edge_set_;
    public:
        BarcodeDistGraphLabeler(const G &g, const BarcodeMapper &mapper, const edge_set_t &edge_set) :
                base(g), barcode_mapper_(mapper), edge_set_(edge_set) { }

        size_t barcode_threshold = 5;

        std::string label(EdgeId e) const {
            std::string ans;
            std::vector <std::string> head_labels;
            std::vector <std::string> tail_labels;
            for (auto edge : edge_set_) {
                if (barcode_mapper_.IntersectionSize(e, edge) >= barcode_threshold) {
                    std::string str = ToString(this->graph().int_id(edge)) + ": " +
                                      std::to_string(barcode_mapper_.IntersectionSize(e, edge)) + ", ";
                    head_labels.push_back(str);
                }
                if (barcode_mapper_.IntersectionSize(edge, e) >= barcode_threshold) {
                    std::string str = ToString(this->graph().int_id(edge)) + ": " +
                                      std::to_string(barcode_mapper_.IntersectionSize(edge, e)) + ", ";
                    tail_labels.push_back(str);
                }
            }
            ans += ("Id: " + std::to_string(this->graph().int_id(e)) + ' ' +  "Length: " + 
                std::to_string(this->graph().length(e)) + ' ' + "Coverage: " + std::to_string(this->graph().coverage(e)) + '\n');
            // INFO(std::to_string(barcode_mapper_.size("head")));
            // ans += (std::to_string(barcode_mapper_.size("head")) + '\n');
            ans += label_from_vector(head_labels, "head");
            ans += label_from_vector(tail_labels, "tail");
            return ans;
        }

        string label_from_vector (const vector <string>& labels, const string& prefix) const {
            std::string ans;
            ans += prefix + ":\n";
            int buffer_len = 0;
            for (auto label : labels) {
                if (buffer_len < 6) {
                    ans += label;
                    buffer_len++;
                }
                else {
                    ans += '\n';
                    buffer_len = 0;
                }
            }
            ans += '\n';
            return ans;
        }
    };

    template <class Graph>
    class TslrVisualizer {
        typedef BarcodeMapper b_mapper;
    private:
        const conj_graph_pack &gp_;
        const b_mapper &barcode_mapper_;
    public:
        TslrVisualizer(const conj_graph_pack &gp, const b_mapper &barcode_mapper) :
                gp_(gp), barcode_mapper_(barcode_mapper) { }

        std::string pics_folder = cfg::get().output_dir + '/' + "pictures";

        void DrawEdgeComponent(const EdgeId &edge, size_t num) {
            size_t max_vertices = 60;
            size_t edge_length_bound = 10000;
            auto component = omnigraph::EdgeNeighborhood(gp_.g, edge, max_vertices, edge_length_bound);
            auto edge_set = component.edges();
            omnigraph::EdgePosGraphLabeler <Graph> pos_labeler(gp_.g, gp_.edge_pos);
            BarcodeDistGraphLabeler<Graph> barcode_labeler(gp_.g, barcode_mapper_, edge_set);
            omnigraph::CompositeLabeler <Graph> composite_labeler(pos_labeler, barcode_labeler);
            auto colorer = omnigraph::visualization::DefaultColorer(gp_.g);
            omnigraph::visualization::WriteComponent(component,
                                                     pics_folder + "/vertex_component_" + std::to_string(num) + "_.dot",
                                                     colorer, composite_labeler);
        }

        void DrawGraph() {
            omnigraph::GraphComponent <Graph> all(gp_.g, gp_.g.begin(), gp_.g.end());
            omnigraph::EdgePosGraphLabeler <Graph> pos_labeler(gp_.g, gp_.edge_pos);
            std::set <EdgeId> edge_set;
            edge_it_helper helper(gp_.g);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                edge_set.insert(*it);
            }
            BarcodeDistGraphLabeler<Graph> barcode_labeler(gp_.g, barcode_mapper_, edge_set);
            omnigraph::CompositeLabeler <Graph> composite_labeler(pos_labeler, barcode_labeler);
            auto colorer = omnigraph::visualization::DefaultColorer(gp_.g);
            omnigraph::visualization::WriteSimpleComponent(all, pics_folder + "/nc_full_graph.dot", colorer, composite_labeler);
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
