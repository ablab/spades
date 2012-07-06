#pragma once

#include "graph_pack.hpp"
#include "omni/visualization_utils.hpp"


namespace debruijn_graph {


class OnlineVisualizer {
private:
	typedef conj_graph_pack graph_pack;
	typedef graph_pack::graph_t Graph;
	typedef Graph::VertexId VertexId;
	typedef Graph::EdgeId EdgeId;

	const conj_graph_pack &gp_;
	size_t picture_counter_;
	string folder_;
	string file_name_base_;
	size_t max_vertices_;
	size_t edge_length_bound_;
	map<EdgeId, string> coloring_;
	EdgesPositionHandler<Graph> positions_;
	total_labeler_graph_struct graph_struct_;
	total_labeler tot_lab_;

	class CountingCallback : public PathProcessor<Graph>::Callback {
	private:
		size_t cnt;
	public:
		CountingCallback() : cnt(0) {
		}

		virtual void HandlePath(const vector<EdgeId>& path) {
			cnt++;
		}

		size_t get_cnt() const {
			return cnt;
		}
	};

public:
	OnlineVisualizer(const conj_graph_pack &gp) : gp_(gp), picture_counter_(0), folder_("pictures"), file_name_base_("picture"),
		max_vertices_(40), edge_length_bound_(1000), positions_(gp.g, 0, true, 15), graph_struct_(gp.g, &gp.int_ids, &positions_), tot_lab_(&graph_struct_) {
		NewPathColorer<Graph> colorer(gp.g);
		NewExtendedSequenceMapper<Graph> mapper(gp_.g, gp_.index, gp_.kmer_mapper, gp_.k_value + 1);
		MappingPath<EdgeId> path1 = mapper.MapSequence(gp.genome);
		MappingPath<EdgeId> path2 = mapper.MapSequence(!gp.genome);
		coloring_ = colorer.ColorPath(path1.simple_path(), path2.simple_path());
		PosFiller<Graph, NewExtendedSequenceMapper<Graph> > filler(gp.g, mapper, positions_);
		filler.Process(gp.genome, "ref0");
		filler.Process(!gp.genome, "ref1");
	}

	void run() {
		while(true) {
			cout << "> ";
			string command;
			getline(cin, command);
			stringstream ss(command);
			string com;
			ss >> com;
			if(com == "exit")
				break;
			else if(com == "set_folder") {
				SetFolder(ss);
			} else if(com == "set_file_name") {
				SetFileName(ss);
			} else if(com == "set_max_vertices") {
				SetMaxVertices(ss);
			} else if(com == "fill_pos") {
				FillPos(ss);
			} else if(com == "clear_pos") {
				ClearPos(ss);
			} else if(com == "vertex") {
				DrawVertex(ss);
			} else if(com == "edge") {
				DrawEdge(ss);
			} else if(com == "position") {
				DrawGenomePosition(ss);
			} else if(com == "paths") {
				FindPaths(ss);
			} else {
				cout << "ignoring command " << command << endl;
			}
		}
	}

	void FindPaths(stringstream & ss) {
		size_t from;
		size_t to;
		size_t max_length;
		ss >> from;
		ss >> to;
		ss >> max_length;
		CountingCallback callback;
		PathProcessor<Graph> pp(gp_.g, 0, max_length, gp_.g.int_ids().ReturnVertexId(from), gp_.g.int_ids().ReturnVertexId(to), callback);
		pp.Process();
		cout << callback.get_cnt() << " paths found" << endl;
	}

	void SetMaxVertices(stringstream & ss) {
		size_t max;
		ss >> max;
		max_vertices_ = max;
	}

	void FillPos(stringstream & ss) {
		string name;
		string file;
		ss >> name;
		ss >> file;
		if(!fileExists(file)) {
			cout << "file " << file << " does not exist" << endl;
			return;
		}
		NewExtendedSequenceMapper<Graph> mapper(gp_.g, gp_.index, gp_.kmer_mapper, gp_.k_value);
		PosFiller<Graph, NewExtendedSequenceMapper<Graph>> filler(gp_.g, mapper, positions_);
		io::Reader irs(file);
		while(!irs.eof()){
			io::SingleRead read;
			irs >> read;
			if (read.IsValid()) {
				Sequence contig = read.sequence();
				filler.Process(contig, name + "_" + read.name());
				filler.Process(!contig, name + "_" + read.name() + "_RC");
			}
		}
	}

	void ClearPos(stringstream & ss) {
		positions_.clear();
		NewExtendedSequenceMapper<Graph> mapper(gp_.g, gp_.index, gp_.kmer_mapper, gp_.k_value);
		PosFiller<Graph, NewExtendedSequenceMapper<Graph>> filler(gp_.g, mapper, positions_);
		filler.Process(gp_.genome, "ref0");
		filler.Process(!gp_.genome, "ref1");
	}

	void SetFolder(stringstream & ss) {
		ss >> folder_;
	}

	void SetFileName(stringstream & ss) {
		ss >> file_name_base_;
	}

	void DrawVertex(stringstream & ss) {
		size_t id;
		ss >> id;
		DrawPicture(gp_.g.int_ids().ReturnVertexId(id));
	}

	void DrawEdge(stringstream & ss) {
		size_t id;
		ss >> id;
		DrawPicture(gp_.g.int_ids().ReturnEdgeId(id));
	}

	void DrawGenomePosition(stringstream & ss) {
		int position;
		ss >> position;
		Sequence all = gp_.genome;
		if(!ss.eof()) {
			string param;
			ss >> param;
			if(param == "--rc") {
				all = !all;
			}
		}
		if(position + gp_.k_value + 1 > all.size())
			cout << "Ignoring request. Position is out of range : required position is " << position << " while length of the sequence is " << all.size() << endl;
		else
			DrawPicture(all.Subseq(position).start<runtime_k::RtSeq::max_size>(gp_.k_value + 1));
	}

	void DrawPicture(runtime_k::RtSeq kmer) {
		kmer = gp_.kmer_mapper.Substitute(kmer);
		if(!gp_.index.contains(kmer)) {
			cout << "No corresponding graph location " << endl;
			return;
		}
		pair<EdgeId, size_t> position = gp_.index.get(kmer);
		if(position.second * 2 < gp_.g.length(position.first))
			DrawPicture(gp_.g.EdgeStart(position.first));
		else
			DrawPicture(gp_.g.EdgeEnd(position.first));
	}

	void DrawPicture(EdgeId e) {
		DrawPicture(gp_.g.EdgeStart(e));
	}

	void DrawPicture(VertexId v) {
		make_dir(folder_);
		stringstream namestream;
		namestream << folder_ << "/" << file_name_base_  << "_" << picture_counter_ << ".dot";
		string file_name = namestream.str();
		stringstream linksteam;
		linksteam << folder_ << "/" << file_name_base_ << "_latest.dot";
		VertexNeighborhoodFinder<Graph> splitter(gp_.g, v, max_vertices_, edge_length_bound_);
		EdgePosGraphLabeler<Graph> labeler(gp_.g, gp_.edge_pos);
		WriteComponents<Graph>(gp_.g, splitter, file_name, *DefaultColorer(gp_.g, coloring_), tot_lab_);
		WriteComponents<Graph>(gp_.g, splitter, linksteam.str(), *DefaultColorer(gp_.g, coloring_), tot_lab_);
		cout << "Picture written to " << file_name << endl;
		picture_counter_++;
	}

};

}
