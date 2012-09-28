#pragma once

namespace cap {

const size_t kDefaultMaxColorsUsed = 100;

typedef size_t TColor;

class TColorSet {
	typedef std::bitset <kDefaultMaxColorsUsed> TBitSet;

private:
	TBitSet bitset_; 

public:
	TColorSet() : bitset_ () {
	}

	TColorSet(const string &bitset_string) : bitset_ (bitset_string) {
	}

	TColorSet(const unsigned long bitset_ul) : bitset_ (bitset_ul) {
	}

	TColorSet(const TBitSet &base_set) {
		bitset_ = base_set;
	}

	const TBitSet &getBitset() const {
		return bitset_;
	}

	inline TColorSet operator | (const TColorSet &other) const {
		return TColorSet(getBitset() | other.getBitset());
	}

	inline bool operator == (const TColorSet &other) const {
		return getBitset() == other.getBitset();
	}

	inline bool operator != (const TColorSet &other) const {
		return getBitset() != other.getBitset();
	}

	bool operator [] (size_t pos) const {
		return bitset_[pos];
	}
	
	bool operator < (const TColorSet &other) const {
		const TBitSet &other_bitset = other.getBitset();
		for (size_t i = 0; i < kDefaultMaxColorsUsed; ++i) {
			if (bitset_[i] != other_bitset[i]) {
				return bitset_[i] < other_bitset[i];
			}
		}
		return false;
	}

	bool any() const {
		return bitset_.any();
	}

	void SetBit(size_t bit_number, bool value) {
		bitset_[bit_number] = value;
	}

	string ToString() {
		return bitset_.to_string();
	}

	static TColorSet SingleColor(TColor color) {
		TBitSet bitset;
		bitset[color] = 1;
		
		return TColorSet(bitset);
	}
};


// ColorGenerator: Singleton class for generating colors
// ColorGenerator generates max_colors different colors in HSV format
// First color is always black
class ColorGenerator {
	 size_t max_colors_;

	// Hue array of needed size
	 vector <double> hue_array_;

	 static double GenerateIthColor(const size_t color_number) {
		double hue_value = 0;
		int accumulated_exp = 0;

		for (size_t i = 0; (1ul << i) <= color_number; ++i) {
			bool bit = (color_number >> i) & 1;
			if (bit) {
				hue_value = hue_value / (1 << accumulated_exp);
				hue_value += 0.5;
				accumulated_exp = 0;
			} else {
				accumulated_exp++;
			}
		}

		return hue_value;
	}

public:
	ColorGenerator(const size_t max_colors = kDefaultMaxColorsUsed) : max_colors_(0), hue_array_() {
		GenerateColors(max_colors);
	}

  void GenerateColors(const size_t number_of_colors) {
		// If all needed colors were already generated, do nothing
		if (number_of_colors <= max_colors_) {
			return;
		}

		hue_array_.resize(number_of_colors);
		for (size_t i = max_colors_; i < number_of_colors; ++i) {
			hue_array_[i] = GenerateIthColor(i);
		}
		max_colors_ = number_of_colors;
	}

	string GetIthColor(const size_t color_number) const {
		VERIFY(color_number < max_colors_);

		// black one is the very special
		if (color_number == 0) {
			return "0 0 0";
		}

		return str(
			boost::format("%.3lf %.3lf %.3lf") % hue_array_[color_number - 1] % 1 % 1
			);
	}

  static ColorGenerator instance() {
    static ColorGenerator instance;
    return instance;
  }

};

template<class Graph, class Element>
class ElementColorHandler: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;

	// For each element will store a bitmask of used there colors.
	restricted::map<Element, TColorSet > data_;

	// Maximum number of different colors that may be used in coloring
	size_t max_colors_;

public:
	// here we have no VERIFYcation. However, there is in color generator.
	 string color_str(const TColor color) const {
		return ColorGenerator::instance().GetIthColor((size_t) color);
	}

	string color_str(const TColorSet &color_set) const {
		if (!color_set.any()) {
			return color_str((TColor) 0);
		}
		string result = "";
		for (size_t i = 0; i < max_colors_; ++i) {
			if (!color_set[i]) continue;
			if (result.length() != 0) {
				result += ':';
			}
			result += color_str((TColor) (i + 1));
		}
		return result;
	}

	ElementColorHandler(const Graph& g, const size_t max_colors = kDefaultMaxColorsUsed) :
    base(g, "ElementColorHandler"),
    max_colors_(max_colors) {
	}

	void PaintElement(Element e, const TColor color) {
		TColorSet &e_colors = data_[e];
		e_colors.SetBit((size_t) color, 1);
	}

	void PaintElement(Element e, const TColorSet &color_set) {
		TColorSet &e_colors = data_[e];
		e_colors = e_colors | color_set;
	}

	TColorSet Color(Element e) const {
		auto it = data_.find(e);
		if (it == data_.end())
			return TColorSet();
		else
			return it->second;
	}

	string ColorStr(Element e) const {
		return color_str(Color(e));
	}

	/*virtual*/
	void HandleDelete(Element e) {
		data_.erase(e);
	}

	void Erase(Element e) {
		data_.erase(e);
	}
};

template<class Graph>
class ColorHandler: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	ElementColorHandler<Graph, EdgeId> edge_color_;
	ElementColorHandler<Graph, VertexId> vertex_color_;
public:

	ColorHandler(const Graph& g, const size_t max_colors = kDefaultMaxColorsUsed) :
			base(g, "ColorHandler"),
			edge_color_(g, max_colors),
			vertex_color_(g, max_colors) {

	}

	void PaintEdge(EdgeId e, const TColor color) {
		edge_color_.PaintElement(e, color);
	}

	void PaintEdge(EdgeId e, const TColorSet &color_set) {
		edge_color_.PaintElement(e, color_set);
	}

	void PaintVertex(VertexId v, const TColor color) {
		vertex_color_.PaintElement(v, color);
	}

	void PaintVertex(VertexId v, const TColorSet &color_set) {
		vertex_color_.PaintElement(v, color_set);
	}

	TColorSet Color(EdgeId e) const {
		return edge_color_.Color(e);
	}

	TColorSet Color(VertexId v) const {
		return vertex_color_.Color(v);
	}

	map<EdgeId, string> EdgeColorMap() const {
		map<EdgeId, string> answer;
		for (auto it = this->g().SmartEdgeBegin(); !it.IsEnd(); ++it) {
			answer[*it] = edge_color_.ColorStr(*it);
		}
		return answer;
	}

	map<VertexId, string> VertexColorMap() const {
		map<VertexId, string> answer;
		for (auto it = this->g().begin(); it != this->g().end(); ++it) {
			answer[*it] = vertex_color_.ColorStr(*it);
		}
		return answer;
	}

	/*virtual*/
	void HandleDelete(EdgeId e) {
		edge_color_.Erase(e);
	}

	/*virtual*/
	void HandleDelete(VertexId v) {
		vertex_color_.Erase(v);
	}

	/*virtual*/
	void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
		VERIFY(old_edges.size() > 0);
//		auto color = Color(old_edges.front());
		for (auto it = old_edges.begin(); it != old_edges.end(); ++it) {
//			VERIFY(color == Color(*it));
			PaintEdge(new_edge, Color(*it));
		}
//		Paint(new_edge, color);
	}

	/*virtual*/
	void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		//todo temporary verification
//		VERIFY(Color(edge2) == edge_type::black && new_edge == edge2);
		PaintEdge(new_edge, Color(edge2));
		PaintEdge(new_edge, Color(edge1));
	}

	/*virtual*/
	void HandleSplit(EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge_2) {
		PaintVertex(this->g().EdgeEnd(new_edge_1), Color(old_edge));
		PaintEdge(new_edge_1, Color(old_edge));
		PaintEdge(new_edge_2, Color(old_edge));
	}
};


template<class Graph>
void SaveColoring(const Graph& g
		, const IdTrackHandler<Graph>& int_ids
		, const ColorHandler<Graph>& coloring
		, const string& filename) {
	GraphComponent<Graph> whole_graph(g);
	ofstream stream((filename + ".clr").c_str());
	stream << whole_graph.v_size() << endl;
	for (auto it = whole_graph.v_begin(); it != whole_graph.v_end(); ++it) {
		stream << g.int_id(*it) << " " << coloring.Color(*it).ToString() << endl;
	}
	stream << whole_graph.e_size() << endl;
	for (auto it = whole_graph.e_begin(); it != whole_graph.e_end(); ++it) {
		stream << g.int_id(*it) << " " << coloring.Color(*it).ToString() << endl;
	}
}

template<class Graph>
void LoadColoring(const Graph& g
		, const IdTrackHandler<Graph>& int_ids
		, ColorHandler<Graph>& coloring
		, const string& filename) {
	ifstream stream((filename + ".clr").c_str());
	size_t v_count;
	stream >> v_count;
	for (size_t i = 0; i < v_count; ++i) {
		size_t id;
		stream >> id;
		string color_string;
		stream >> color_string;
		coloring.Paint(int_ids.ReturnVertexId(id), TColorSet(color_string));
	}
	size_t e_count;
	stream >> e_count;
	for (size_t i = 0; i < e_count; ++i) {
		size_t id;
		stream >> id;
		string color_string;
		stream >> color_string;
		coloring.Paint(int_ids.ReturnEdgeId(id), TColorSet(color_string));
	}
}

template<class Graph>
auto_ptr<GraphColorer<Graph>> ConstructColorer(
		const ColorHandler<Graph>& coloring) {
	return auto_ptr<GraphColorer<Graph>>(
			new CompositeGraphColorer<Graph>(
					new MapColorer<typename Graph::VertexId>(
							coloring.VertexColorMap()),
					new MapColorer<typename Graph::EdgeId>(
							coloring.EdgeColorMap())));
}

template<class Graph>
auto_ptr<GraphColorer<Graph>> ConstructBorderColorer(const Graph& g,
		const ColorHandler<Graph>& coloring) {
	return auto_ptr<GraphColorer<Graph>>(
			new CompositeGraphColorer<Graph>(
					new BorderVertexColorer<Graph>(g),
					new MapColorer<typename Graph::EdgeId>(
							coloring.EdgeColorMap())));
}

// Temporary while have only two colors
TColor kRedColor = (TColor) 0;
TColor kBlueColor = (TColor) 1;
TColorSet kRedColorSet = TColorSet::SingleColor(kRedColor);
TColorSet kBlueColorSet = TColorSet::SingleColor(kBlueColor);
TColorSet kVioletColorSet = kRedColorSet | kBlueColorSet;

}

