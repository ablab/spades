#ifndef REPEAT_HPP
#define REPET_HPP

namespace debruijn_graph {

	typedef enum { TOPOLOGY, LENGTH, PAIREDINFO } kind_of_repeat;

	bool if_unique_val(double val, const vector<double>& ar, double epsilon=0.000001) {
		unsigned counter = 0;
		for (auto it = ar.begin(); it != ar.end(); ++it) {
			if (fabs(val - *it)<=epsilon) counter+=1;
		}
		return (counter==1);
	}

	double get_max(const vector<vector<double>>& vec, int& max_id_i, int& max_id_j, double epsilon=0.000001) {
		double max_val = 0;
		max_id_i = -1;
		max_id_j = -1;
		for (unsigned i = 0; i < vec.size(); ++i) {
			for (unsigned j = 0; j < vec[i].size(); ++j) {
				if (vec[i][j] - max_val > epsilon ){
					max_val = vec[i][j];
					max_id_i = i;
					max_id_j = j;
				}
			}
		}
		return max_val;
	}
	
	void hide(vector<vector<double>>& vec, int i, double default_val = -1) {
		for (unsigned j = 0; j < vec[i].size(); ++j) vec[i][j] = default_val;
	}

	void hide(vector<vector<double>>& vec, int i, int j, double default_val = -1) {
		for (unsigned k = 0; k < vec[i].size(); ++k) vec[i][k] = default_val;
		for (unsigned k = 0; k < vec.size(); ++k) vec[k][j] = default_val;
	}

	template <class graph_pack>
	class Repeat {
		
		typedef typename Graph::EdgeId EdgeId;
		typedef typename Graph::VertexId VertexId;

		const graph_pack &gp_;
		vector<EdgeId> incoming_edges_;
		vector<EdgeId> outgoing_edges_;
		vector<EdgeId> component_;
		const double repeat_length_upper_threshold_;
		const map<EdgeId, kind_of_repeat> edge_to_kind_;
		FILE* file;



		template <class EdgesPositionHandlerT> 
		bool match( const vector<EdgeId> &path, const unsigned current_id, const unsigned current_start,
			EdgesPositionHandlerT& ref_pos) {
			if (path.size() == current_id ){
				return true;
			}
			EdgeId edge = path[current_id];
			auto pos_it = ref_pos.edges_positions().find(edge);
			if ( current_id == path.size() - 1 && pos_it->second.size() > 1 ){
				return false;
			}
			bool matched = false;
			for (size_t i = 0; i < pos_it->second.size(); ++i) {
				auto start = pos_it->second[i].start();
				if ( fabs(start - current_start) < 2 ) {
					auto end = pos_it->second[i].end();
					matched = match( path, current_id + 1, end + 1, ref_pos);
				}
			}
			return matched;
		}

		bool MatchReference( const vector<EdgeId>& path) {
			auto ref_pos = gp_.edge_pos;
			EdgeId edge = path[0];
			auto pos_it = ref_pos.edges_positions().find(edge);
			if ( pos_it->second.size() == 1 ){
				auto next_start = pos_it->second[0].end();
				return match( path, 1, next_start + 1, ref_pos );
			}
			return false;
		}

		bool IfContainsLinearlyDependentRows(const vector<vector<double>>& transition_probabilities ) const {
			for (unsigned i = 0; i < transition_probabilities.size() - 1; ++i) {
				bool dependent = true;
				for (unsigned k = i+1; k < transition_probabilities.size(); ++k) 
					for (unsigned j = 0; j < transition_probabilities.size(); ++j) {
						if (fabs(transition_probabilities[i][j] - transition_probabilities[k][j]) > 0.001) { 
							dependent = false;
							break;
						}
					}
				if (dependent) return true;
			}
			return false;
		}
	
		void ChoosePairsGreedy(vector<vector<double>>& transition_probabilities, vector<pair<EdgeId,EdgeId>>& pairs_of_edges, double match_quality_threshold) {
			if (IfContainsLinearlyDependentRows (transition_probabilities)) 
				return;
			unsigned counter(0);
			while (counter < transition_probabilities.size()) {
				int max_id_i(0), max_id_j(0);
				double max_val = get_max(transition_probabilities, max_id_i, max_id_j);
				if (max_id_i != -1) {
					if (!if_unique_val(max_val, transition_probabilities[max_id_i])) {
						hide(transition_probabilities, max_id_i);
					}
					else {
						hide(transition_probabilities, max_id_i, max_id_j);
						if (max_val > match_quality_threshold) {
							pairs_of_edges.push_back(make_pair(incoming_edges_[max_id_i], outgoing_edges_[max_id_j]));
							fprintf(file,"pair: %d %d %5.4f\n", max_id_i, max_id_j, max_val);
						}
					}
				}
				counter += 1;
			}
		}

		//TODO: Dima: make gp_.edge_pos.IsConsistentWithGenome() a const method
		void ChooseConsistentPairs( const vector< vector <double> >& transition_probabilities ) {
			int i = 0, j = 0;
			for (auto in_edge = incoming_edges_.begin(); in_edge != incoming_edges_.end(); ++in_edge, ++i) {
				for ( auto out_edge = outgoing_edges_.begin(); out_edge != outgoing_edges_.end(); ++out_edge, ++j ) {
					vector<EdgeId> path;
					path.push_back(*in_edge);
					path.insert(path.begin()+1, component_.begin(), component_.end());
					path.push_back(*out_edge);
					DEBUG("check if consistent with genome..");
					/*if (gp_.edge_pos.IsConsistentWithGenome(path)) {
						DEBUG("pair: " << transition_probabilities[i][j] << "\n");
					}*/
				}
			}
		}

		void bfs ( const EdgeId& edge,  set<EdgeId>& visited_edges, int& curLen, int& maxPathLen) {
			visited_edges.insert(edge);
			auto incoming_edges = gp_.g.IncomingEdges(gp_.g.EdgeStart(edge));
			for ( auto e = incoming_edges.begin(); e != incoming_edges.end(); ++e) {
				if ( find(component_.begin(), component_.end(), *e) != component_.end() && visited_edges.find(*e) == visited_edges.end() ){
					curLen += gp_.g.length(*e);
					if (curLen > maxPathLen) maxPathLen = curLen;
					bfs(*e, visited_edges, curLen, maxPathLen);
				}
			}
			auto outgoing_edges_ = gp_.g.OutgoingEdges(gp_.g.EdgeEnd(edge));
			for ( auto e = outgoing_edges_.begin(); e != outgoing_edges_.end(); ++e) {
				if ( find(component_.begin(), component_.end(), *e) != component_.end() && visited_edges.find(*e) == visited_edges.end() ){
					curLen += gp_.g.length(*e);
					if (curLen > maxPathLen) maxPathLen = curLen;
					bfs(*e, visited_edges, curLen, maxPathLen);
				}
			}
		}		
	
		public:
		explicit Repeat(const graph_pack& gp, const vector<EdgeId>& incoming_edges, const vector<EdgeId>& outgoing_edges, const vector<EdgeId>& component, double repeat_length_upper_threshold,
				const map<EdgeId, kind_of_repeat>& edge_to_kind, FILE* out_file) : gp_(gp),
										incoming_edges_(incoming_edges),
										outgoing_edges_(outgoing_edges),
										component_(component), 
										repeat_length_upper_threshold_(repeat_length_upper_threshold),
										edge_to_kind_(edge_to_kind),
										file(out_file) {}

			      
		template<class EdgeQualityLabeler>
		bool IfContainsOnlyGenomicEdges( const EdgeQualityLabeler& quality_labeler ) const {
			for (auto iter = component_.begin(); iter != component_.end(); ++iter) {
		;		if (quality_labeler.quality(*iter) < 0.5) {
					return false;
				}
			}		
			for ( auto e = incoming_edges_.begin(); e != incoming_edges_.end(); ++e ) {
				if (quality_labeler.quality(*e) < 0.5) {
					return false;
				}
			}
			for ( auto e = outgoing_edges_.begin(); e != outgoing_edges_.end(); ++e ) {
				if (quality_labeler.quality(*e) < 0.5) {
					return false;
				}
			}
			return true;
		}
			
		template<class EdgeQualityLabeler>
		bool IfRepeatByQuality( const EdgeQualityLabeler& quality_labeler ) const {

			for ( auto e = component_.begin(); e != component_.end(); ++e ) {
				if ( quality_labeler.quality(*e) <= 1.5 ) return false;
			}
			for ( auto e = incoming_edges_.begin(); e != incoming_edges_.end(); ++e ) {
				 if ( quality_labeler.quality(*e) > 1.5 ) return false;
			}
			for ( auto e = outgoing_edges_.begin(); e != outgoing_edges_.end(); ++e ) {
				 if ( quality_labeler.quality(*e) > 1.5 ) return false;
			}
			return true;
		}

		template <class DetailedCoverage, class EdgeQualityLabeler>
		void GetComponentInfo(const DetailedCoverage& coverage, const vector< vector <double> >& transition_probabilities, const EdgeQualityLabeler& quality_labeler ) const  {
			fprintf(file,"Component: \n");
			for ( auto iter = component_.begin(); iter != component_.end(); ++iter ) {
				fprintf(file,"%lu edge length: %lu average edge coverage %5.4f quality: %5.2f",gp_.g.int_id(*iter), gp_.g.length(*iter), gp_.g.coverage(*iter),
					quality_labeler.quality(*iter));
				auto repeat_type = edge_to_kind_.find(*iter);
				VERIFY(repeat_type != edge_to_kind_.end()); 
				if ( repeat_type->second == TOPOLOGY){
					fprintf(file,"TOPOLOGY\n");
				}
				else if  (repeat_type->second == LENGTH ){
					fprintf(file,"LENGTH\n");
				}
				else if (repeat_type->second == PAIREDINFO ){
					fprintf(file,"PAIREDINFO\n");
				}
			 }
			fprintf(file, "incoming edges:\n");
			for ( auto iter = incoming_edges_.begin(); iter != incoming_edges_.end(); ++iter ) {
			 	fprintf(file, "%lu edge length: %lu outgoing edge coverage: %5.4f average edge coverage %5.4f quality: %5.2f\n", gp_.g.int_id(*iter), gp_.g.length(*iter), 
					coverage.GetOutCov(*iter), gp_.g.coverage(*iter), quality_labeler.quality(*iter));
			}
			fprintf(file,"outgoing edges: \n");
			for ( auto iter = outgoing_edges_.begin(); iter != outgoing_edges_.end(); ++iter ) {
			 	fprintf(file, "%lu edge length: %lu incoming edge coverage: %5.4f average edge coverage %5.4f quality: %5.2f\n", gp_.g.int_id(*iter), gp_.g.length(*iter), 
					coverage.GetInCov(*iter), gp_.g.coverage(*iter), quality_labeler.quality(*iter));
			}
			bool correct_component = IfRepeatByQuality( quality_labeler  );
			if (!correct_component) {
				fprintf(file,"repeat is detected incorrectly\n");
			}
			if (transition_probabilities.size() > 0) {
			 	fprintf(file,"transition probabilities:\n");
				for (auto vec = transition_probabilities.begin(); vec != transition_probabilities.end(); ++vec ) {
					for ( auto prob = vec->begin(); prob != vec->end(); ++prob ) {
						fprintf(file,"%5.4f ",*prob);
					}
					fprintf(file,"\n");
				}
			}
				
		}

		bool IfContainsOnlyShortEdges() const {
			for ( auto it = component_.begin(); it != component_.end(); ++it ) {
				auto repeat_type = edge_to_kind_.find(*it);
				VERIFY(repeat_type != edge_to_kind_.end()); 
				if (gp_.g.length(*it) >= repeat_length_upper_threshold_ || repeat_type->second == TOPOLOGY )
					return false;
			}
			return true;
		}

		void CountDistance( const EdgeId& edge_in, const EdgeId& edge_out, int& distance ) const {
			// gets a repetitive component and calculates the length of the longest path in it
			if ( gp_.g.EdgeEnd(edge_in) == gp_.g.EdgeStart(edge_out) ) return;
			for ( auto edge = component_.begin(); edge != component_.end(); ++edge ){
				if ( gp_.g.EdgeEnd(edge_in) == gp_.g.EdgeStart(*edge) ) {
					distance += gp_.g.length(*edge);
					CountDistance( *edge, edge_out, distance);
					return;
				}
			}
		}

		// the length of the longest path in kmers
		int GetLongestPathLength() const {
			set<EdgeId> visited_edges;
			vector<vector<EdgeId>> paths;
			int maxPathLen = gp_.g.length(component_[0]);
			for ( auto edge = component_.begin(); edge != component_.end(); ++edge ){
				if (visited_edges.find(*edge) != visited_edges.end()) continue;
				int curLen = gp_.g.length(*edge);
				visited_edges.insert(*edge);
				 bfs(*edge, visited_edges, curLen, maxPathLen);
			}
			return maxPathLen;
		}

		void dfs( const VertexId& vStart, const VertexId& vEnd, set<EdgeId>& visited,  vector<EdgeId>& path) const {
			for ( auto edge = component_.begin(); edge != component_.end(); ++edge ){
				if (visited.find(*edge) != visited.end()){
					continue;
				}
				if ( vStart == gp_.g.EdgeStart(*edge) ){
					if ( vEnd == gp_.g.EdgeEnd(*edge) ){
						 path.push_back(*edge);
						 visited.insert(*edge);
						 return;
					}
					visited.insert(*edge);
					dfs( gp_.g.EdgeEnd(*edge), vEnd, visited, path );
					path.insert(path.begin(), *edge);
					return;
				}
			}
			return;
		}

		void SetPaths(const vector<pair<EdgeId,EdgeId>>& pairs_of_edges, vector<vector<EdgeId>>& resolved_paths) {
			for (unsigned i = 0; i < pairs_of_edges.size(); ++i) {
				EdgeId in_edge = pairs_of_edges[i].first;
				EdgeId out_edge = pairs_of_edges[i].second;
				vector<EdgeId> path;
				set<EdgeId> visited;
				dfs(gp_.g.EdgeEnd(in_edge), gp_.g.EdgeStart(out_edge), visited, path);
				path.insert(path.begin(),in_edge);
				path.push_back(out_edge);
				if (MatchReference(path)) { fprintf(file,"match!\n"); }
				else {fprintf(file,"does not match!\n");}
				resolved_paths.push_back(path);
			}

		}

		template <class DetailedCoverage, class EdgeQualityLabeler, class KmerIndex>
		bool Resolve( BucketMapper<Graph, KmerIndex> &bm, const DetailedCoverage& coverage,
							const EdgeQualityLabeler& quality_labeler,
							vector<vector<EdgeId>>& resolved_paths) { 
			vector<pair<EdgeId,EdgeId>> pairs_of_edges;
			vector< vector <double> > transition_probabilities ;
			for ( unsigned i = 0; i < incoming_edges_.size(); i++) {
				transition_probabilities.push_back(vector<double>(outgoing_edges_.size(),1));
			}
			int in_edge_counter(0);
			for ( auto in_edge = incoming_edges_.begin(); in_edge != incoming_edges_.end(); ++in_edge, ++in_edge_counter ) {
				double in_cov = coverage.GetOutCov(*in_edge);
				int in_bucket = bm.GetCoverageBucket(in_cov);
				int out_edge_counter(0);
				for ( auto out_edge = outgoing_edges_.begin(); out_edge != outgoing_edges_.end(); ++out_edge, ++out_edge_counter ) {
					double out_cov = coverage.GetInCov(*out_edge);
					int out_bucket = bm.GetCoverageBucket(out_cov);
					int distance(0);
					CountDistance(*in_edge, *out_edge, distance);
					double shift = 25;
					double probability = bm.GetProbabilityFromBucketToBucketForDistance (in_bucket, out_bucket, distance, shift) ;
					transition_probabilities[in_edge_counter][out_edge_counter] = probability;
				}
			}
			GetComponentInfo(coverage, transition_probabilities, quality_labeler );
			double match_quality_threshold = 0.1;
			ChoosePairsGreedy(transition_probabilities, pairs_of_edges, match_quality_threshold);
			if (pairs_of_edges.size() > 0) {
				SetPaths(pairs_of_edges, resolved_paths);
				return true;
			}
			return false;
		}
	};

}

#endif
