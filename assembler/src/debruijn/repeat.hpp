#ifndef REPEAT_HPP
#define REPET_HPP

namespace debruijn_graph {

	typedef enum { TOPOLOGY, LENGTH, PAIREDINFO } kind_of_repeat;

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

		void ChoosePairsGreedy( const vector< vector <double> > transition_probabilities, vector<pair<EdgeId,EdgeId>>& pairs_of_edges) {
			unsigned k = 0;			
			set<unsigned> matched;
			DEBUG("MATCH: ");
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
			DEBUG("Component: ");
			for ( auto iter = component_.begin(); iter != component_.end(); ++iter ) {
				DEBUG( gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) <<  " average edge coverage " 
					<< gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << " ");
				auto repeat_type = edge_to_kind_.find(*iter);
				VERIFY(repeat_type != edge_to_kind_.end()); 
				if ( repeat_type->second == TOPOLOGY){
					DEBUG("TOPOLOGY");
				}
				else if  (repeat_type->second == LENGTH ){
					DEBUG("LENGTH");
				}
				else if (repeat_type->second == PAIREDINFO ){
					DEBUG("PAIREDINFO");
				}
			 }
			DEBUG("\n");
			DEBUG("incoming edges: ");
			for ( auto iter = incoming_edges_.begin(); iter != incoming_edges_.end(); ++iter ) {
			 	DEBUG(gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) << " outgoing edge coverage: " << coverage.GetOutCov(*iter) << 
					" average edge coverage " << gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << "\n");
			}
			DEBUG("\n");
			DEBUG("outgoing edges: ");
			for ( auto iter = outgoing_edges_.begin(); iter != outgoing_edges_.end(); ++iter ) {
				DEBUG(gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) << " incoming edge coverage: " << coverage.GetInCov(*iter) <<
					" average edge coverage " << gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << "\n");
			}
			bool correct_component = IfRepeatByQuality( quality_labeler  );
			if (!correct_component) {
				DEBUG("repeat is detected incorrectly\n");
			}
			if (transition_probabilities.size() > 0) {
			 	DEBUG("transition probabilities:\n");
				for (auto vec = transition_probabilities.begin(); vec != transition_probabilities.end(); ++vec ) {
					for ( auto prob = vec->begin(); prob != vec->end(); ++prob ) {
						DEBUG(*prob << " ");
					}
					DEBUG("\n");
				}
			}
			DEBUG("\n");
				
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

		template <class DetailedCoverage>
		bool Resolve( BucketMapper<Graph> &bm, const DetailedCoverage& coverage,
							const EdgeQuality<typename graph_pack::graph_t>& quality_labeler,
							set<vector<EdgeId>>& resolved_paths) { 
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
			if (pairs_of_edges.size() > 0) return true;

			return false;
		}
	};

}

#endif
