#pragma once

namespace cap {

/*
 *  This filter restricts also components where all cops are borschs
 *  (all edges are colored in all colors)
 */
template<class Graph>
class ComponentSingleColorFilter: public GraphComponentFilter<Graph> {
private:
	typedef GraphComponentFilter<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

    ColorHandler<Graph> color_handler_;
    TColorSet restricted_color_;
	size_t max_length_;
	size_t vertex_number_;

public:
	ComponentSingleColorFilter(const Graph &graph, const ColorHandler<Graph> &color_handler,
            const TColorSet &restricted_color, size_t max_length, size_t vertex_number)
        : base(graph),
          color_handler_(color_handler),
          restricted_color_(restricted_color),
          max_length_(max_length),
          vertex_number_(vertex_number) {
	}

	/*virtual*/
	bool Check(const vector<VertexId> &vertices) const {
        TRACE("Check component");
		if (vertices.size() <= vertex_number_)
			return false;

        bool length_flag = false,
             color_flag = false;
		set < VertexId > component(vertices.begin(), vertices.end());
		for (auto iterator = vertices.begin(); iterator != vertices.end();
				++iterator) {
			vector < EdgeId > edges = this->graph().OutgoingEdges(*iterator);
			for (auto edge_iterator = edges.begin();
					edge_iterator != edges.end(); edge_iterator++) {
				if (component.count(this->graph().EdgeEnd(*edge_iterator)) == 1) {
                    if (this->graph().length(*edge_iterator) <= max_length_) {
                        length_flag = true;
                    }
                    if (color_handler_.Color(*edge_iterator) != restricted_color_) {
                        TRACE("Found good color " << color_handler_.Color(*edge_iterator).ToString());
                        color_flag = true;
                    }
                    if (length_flag && color_flag) {
                        return true;
                    }
				}
			}
		}
		return false;
	}
};

template<class Graph>
void PrintColoredGraph(const Graph& g, const ColorHandler<Graph>& coloring,
		const EdgesPositionHandler<Graph>& pos, const string& output_filename) {
	ReliableSplitter<Graph> splitter(g, 30, 1000000);
	LengthIdGraphLabeler<Graph> basic_labeler(g);
	EdgePosGraphLabeler<Graph> pos_labeler(g, pos);

	CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
	WriteComponents(g, splitter, output_filename,
//				*ConstructColorer(coloring),
			*ConstructBorderColorer(g, coloring), labeler);
}

template<class Graph>
void PrintColoredGraphWithColorFilter(const Graph &g, const ColorHandler<Graph> &coloring,
    const EdgesPositionHandler<Graph> &pos, const string &output_filename) {
    
  size_t edge_length_bound = 1000000;
  size_t colors_number = coloring.max_colors();
  TColorSet restricted_color = TColorSet::AllColorsSet(colors_number);

	ReliableSplitter<Graph> splitter(g, 30, edge_length_bound);
	ComponentSingleColorFilter<Graph> filter(g, coloring, restricted_color, edge_length_bound, 2);
	LengthIdGraphLabeler<Graph> basic_labeler(g);
	EdgePosGraphLabeler<Graph> pos_labeler(g, pos);

	CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
	WriteComponents(g, splitter, filter, output_filename,
			*ConstructBorderColorer(g, coloring), labeler);
}

template<class gp_t>
void PrintColoredGraphAlongRef(const gp_t& gp,
		const ColorHandler<Graph>& coloring,
		const string& output_filename) {
	LengthIdGraphLabeler < Graph > basic_labeler(gp.g);
	EdgePosGraphLabeler < Graph > pos_labeler(gp.g, gp.edge_pos);

	CompositeLabeler < Graph > labeler(basic_labeler, pos_labeler);

//		only breakpoints
	TrivialBreakpointFinder<Graph> bp_f(gp.g, coloring, gp.edge_pos);

	WriteComponentsAlongPath(gp.g, bp_f, labeler, output_filename, 1000000,
			30/*000*/, MapperInstance(gp)->MapSequence(gp.genome),
			*ConstructBorderColorer(gp.g, coloring)
//				*ConstructColorer(coloring)
					);
}


}
