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

	
		void ChoosePairsGreedy(vector<vector<double>>& transition_probabilities, vector<pair<EdgeId,EdgeId>>& pairs_of_edges, double match_quality_threshold) {
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
							INFO("pair: " << max_id_i << " " << max_id_j << " " << max_val);
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
				const map<EdgeId, kind_of_repeat>& edge_to_kind) : gp_(gp),
										incoming_edges_(incoming_edges),
										outgoing_edges_(outgoing_edges),
										component_(component), 
										repeat_length_upper_threshold_(repeat_length_upper_threshold),
										edge_to_kind_(edge_to_kind) {}

			        

		bool IfContainsOnlyGenomicEdges( const EdgeQuality<typename graph_pack::graph_t>& quality_labeler ) const {
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
			
		bool IfRepeatByQuality( const EdgeQuality<typename graph_pack::graph_t>& quality_labeler ) const {

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

		template <class DetailedCoverage>
		void GetComponentInfo(const DetailedCoverage& coverage, const vector< vector <double> >& transition_probabilities, const EdgeQuality<typename graph_pack::graph_t>& quality_labeler ) const  {
			INFO("Component: ");
			for ( auto iter = component_.begin(); iter != component_.end(); ++iter ) {
				INFO( gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) <<  " average edge coverage " 
					<< gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << " ");
				auto repeat_type = edge_to_kind_.find(*iter);
				VERIFY(repeat_type != edge_to_kind_.end()); 
				if ( repeat_type->second == TOPOLOGY){
					INFO("TOPOLOGY");
				}
				else if  (repeat_type->second == LENGTH ){
					INFO("LENGTH");
				}
				else if (repeat_type->second == PAIREDINFO ){
					INFO("PAIREDINFO");
				}
			 }
			INFO("\n");
			INFO("incoming edges: ");
			for ( auto iter = incoming_edges_.begin(); iter != incoming_edges_.end(); ++iter ) {
			 	INFO(gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) << " outgoing edge coverage: " << coverage.GetOutCov(*iter) << 
					" average edge coverage " << gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter));
			}
			INFO("\n");
			INFO("outgoing edges: ");
			for ( auto iter = outgoing_edges_.begin(); iter != outgoing_edges_.end(); ++iter ) {
				INFO(gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) << " incoming edge coverage: " << coverage.GetInCov(*iter) <<
					" average edge coverage " << gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter));
			}
			bool correct_component = IfRepeatByQuality( quality_labeler  );
			if (!correct_component) {
				INFO("repeat is detected incorrectly");
			}
			if (transition_probabilities.size() > 0) {
			 	INFO("transition probabilities:\n");
				for (auto vec = transition_probabilities.begin(); vec != transition_probabilities.end(); ++vec ) {
					for ( auto prob = vec->begin(); prob != vec->end(); ++prob ) {
						INFO(*prob << " ");
					}
					INFO("\n");
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
				resolved_paths.push_back(path);
			}

		}

		template <class DetailedCoverage>
		bool Resolve( BucketMapper<Graph> &bm, const DetailedCoverage& coverage,
							const EdgeQuality<typename graph_pack::graph_t>& quality_labeler,
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
			double match_quality_threshold = 0.05;
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
