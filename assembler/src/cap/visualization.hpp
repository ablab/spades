#pragma once

namespace cap {

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
