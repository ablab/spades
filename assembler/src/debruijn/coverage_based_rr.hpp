#ifndef COVERAGE_BASED_RR
#define COVERAGE_BASED_RR

#include <vector>
#include <map>
#include <set>
#include <list>
#include "graph_print_utils.hpp"
#include "indices/perfect_hash_map.hpp"
#include "graph_pack.hpp"
#include "pair_info_improver.hpp"
#include "path_extend/pe_io.hpp"
#include "path_extend/bidirectional_path.hpp"
#include "graphio.hpp"
#include "long_read_storage.hpp"
#include "genome_consistance_checker.hpp"

namespace debruijn_graph{

template <class graph_pack, class EdgeQualityLabeler, class KmerIndex>
class CoverageBasedResolution {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef double coverage_value;

	const graph_pack& gp_;

	const KmerIndex& kmer_index_;
	const EdgeQualityLabeler& quality_labeler_;

	map<EdgeId, kind_of_repeat> edge_to_kind_;
	const double tandem_lower_threshold_;
	const double tandem_upper_threshold_;
	const double repeat_length_upper_threshold_;

	public:
	CoverageBasedResolution( const graph_pack& gpack_arg, const KmerIndex& kmer_index, const EdgeQualityLabeler& quality_labeler,
				double tandem_lower_threshold, double tandem_upper_threshold, double repeat_length_upper_threshold) :

											gp_(gpack_arg),
											kmer_index_(kmer_index),
											quality_labeler_(quality_labeler),
											tandem_lower_threshold_(tandem_lower_threshold), 
											tandem_upper_threshold_(tandem_upper_threshold),
											repeat_length_upper_threshold_(repeat_length_upper_threshold){
	}
	
	
	private:


	bool IfRepeatByPairedInfo( const EdgeId& edge, PairedInfoIndexT<Graph>& clustered_index ) const {
		io::SequencingLibrary<debruijn_config::DataSetData> lib;
		auto improver = PairInfoImprover<Graph>(gp_.g, clustered_index,lib);
		auto inner_map = clustered_index.GetEdgeInfo(edge, 0);
		for (auto I_1 = inner_map.begin(), E = inner_map.end(); I_1 != E; ++I_1) {
			for (auto I_2 = inner_map.begin(); I_2 != E; ++I_2) {
				if (I_1 == I_2) continue;
				EdgeId e1 = (*I_1).first;
				const Point& p1 =
				        *(*I_1).second.begin();//TODO: was const Point& p1 = (*I_1).second, change it
				EdgeId e2 = (*I_2).first;
				const Point& p2 = *(*I_2).second.begin();
				if ( p1.d * p2.d < 0 || p2.d > p1.d ) continue;
				if (!improver.IsConsistent(edge, e1, e2, p1, p2)) {
					return true;
				}
			}
		}
		return false;
	}

	void JoinPaths( const vector<vector<EdgeId> >& paths, 
			const vector<vector<EdgeId>>& resolved_loops,
			vector< vector<EdgeId> >& all_paths ) const {
		map< EdgeId, vector<EdgeId> > start_edge_to_path;
		map< EdgeId, vector<EdgeId> > back_edge_to_path;
		vector< vector<EdgeId> > both_paths;
		for (auto path = resolved_loops.begin(); path != resolved_loops.end(); ++path) {
			both_paths.push_back(*path);
		}
		for (auto path = paths.begin(); path != paths.end(); ++path) {
			both_paths.push_back(*path);
		}
		for (auto path = both_paths.begin(); path != both_paths.end(); ++path) {
			//VERIFY(start_edge_to_path.find(path->front()) == start_edge_to_path.end());
			start_edge_to_path[path->front()] = *path;
		}

		for (auto path = both_paths.begin(); path != both_paths.end(); ++path) {
			//VERIFY(back_edge_to_path.find(path->back()) == back_edge_to_path.end());
			back_edge_to_path[path->back()] = *path;
		}

		for (auto path = both_paths.begin(); path != both_paths.end(); ++path) {
			if ( back_edge_to_path.find( path->front() ) != back_edge_to_path.end()) {
				continue;
			}
			bool updated = true;
			while (updated) {
				auto found_path = start_edge_to_path.find( path->back() );
				if (found_path != start_edge_to_path.end() ) {
					path->insert(path->end(), found_path->second.begin() + 1, found_path->second.end());
				} else {
					updated = false;
				}
			}
			all_paths.push_back(*path);
		}
	}

/*	void dfs( const VertexId& vStart, const VertexId& vEnd, const vector<EdgeId>& component,
		set<EdgeId>& visited,  path_extend::BidirectionalPath& path) const {
		for ( auto edge = component.begin(); edge != component.end(); ++edge ){
			if (visited.find(*edge) != visited.end()){
				continue;
			}
			if ( vStart == gp_.g.EdgeStart(*edge) ){
				if ( vEnd == gp_.g.EdgeEnd(*edge) ){
					path.PushBack(*edge);
					visited.insert(*edge);
					return;
				}
				visited.insert(*edge);
				dfs( gp_.g.EdgeEnd(*edge), vEnd, component, visited, path );
				path.PushFront(*edge);
				return;
			}
		}
		return;
	}
*/
	void visit( const EdgeId& edge, set<EdgeId>& visited_edges, set<VertexId>& grey_vertices,
		vector<EdgeId>& path, const vector<EdgeId>& components, 
		const vector<EdgeId>& singles, vector<EdgeId>& incoming_edges, vector<EdgeId>& outgoing_edges, bool& if_loop ) const {
		VertexId edgeStartVertex = gp_.g.EdgeStart(edge);
		VertexId edgeEndVertex = gp_.g.EdgeEnd(edge);
		if (visited_edges.find(edge) != visited_edges.end() )
			return;
		if ( grey_vertices.find(edgeStartVertex) != grey_vertices.end() ) { 
			if_loop = true;
			return;
		}
		grey_vertices.insert(edgeStartVertex);
		path.push_back(edge);
		visited_edges.insert(edge);
		for ( auto single_edge = singles.begin(); single_edge != singles.end(); ++single_edge ){
			VertexId singleStartVertex = gp_.g.EdgeStart(*single_edge);
			VertexId singleEndVertex = gp_.g.EdgeEnd(*single_edge);
			if ( singleStartVertex == edgeEndVertex ) {
				outgoing_edges.push_back(*single_edge);
			}
			if ( singleEndVertex == edgeStartVertex ) {
				incoming_edges.push_back(*single_edge);
			}
		}
		for ( auto component_edge = components.begin(); component_edge != components.end(); ++component_edge ){
			if ( *component_edge == edge ){
				continue;
			}
			VertexId componentEdgeStartVertex = gp_.g.EdgeStart(*component_edge);
			VertexId componentEdgeEndVertex = gp_.g.EdgeEnd(*component_edge);
			if ( componentEdgeStartVertex == edgeEndVertex || componentEdgeEndVertex == edgeStartVertex ) {
				visit(*component_edge, visited_edges, grey_vertices, path, components, singles, incoming_edges, outgoing_edges, if_loop);
			}
		}
		grey_vertices.insert(edgeStartVertex);
	}


	void GetOtherEdges(path_extend::PathContainer& paths, const set<EdgeId>& used_edges) const {
	// adds edges from the rest of the graph (which was n)
		set<EdgeId> included;
		for (auto iter = gp_.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (used_edges.find(*iter) == used_edges.end() && included.find(*iter) == included.end()){
				paths.AddPair(new path_extend::BidirectionalPath(gp_.g, *iter), new path_extend::BidirectionalPath(gp_.g, gp_.g.conjugate(*iter)));
				included.insert(*iter);
				included.insert(gp_.g.conjugate(*iter));
			}
		}

	}
	void WriteResolved( const set<EdgeId>& used_edges, path_extend::PathContainer& resolved_paths, const string &file_name  ) const {
		path_extend::ContigWriter cw( gp_.g);
		GetOtherEdges( resolved_paths, used_edges );
		cw.writePaths( resolved_paths, file_name );
	}

// class SingleRepeat
//	vector<EdgeId> incoming;
//	vector<EdgeId> outgoing;
//	vector<EdgeId> component;
//	const graph_pack &gp_;
//	bool VerifyComponentByReference;
// 	set<vector<EdgeId> > ResolveRepeat(?);

//	bool VerifyComponentWithUniqueShortEdge (VerifyComponent)?
//  PatchComponentWithUniqueShortEdge (?)

//TODO: excess white lines
//TODO: google styleguide
//TODO setup vim or eclipse via ssh:)
	bool VerifyComponent( vector<EdgeId>& incoming_edges,
			 	vector<EdgeId>& outgoing_edges,
		 		vector<EdgeId>& component ){
		if ( incoming_edges.size() == outgoing_edges.size() ) {
			return true;
		}
		int diff = incoming_edges.size() - outgoing_edges.size();
		if ( diff < 0 ) {
			int counter = 0;
			for ( auto edge = incoming_edges.begin(); edge != incoming_edges.end(); ++edge ) {
				if ( gp_.g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
					counter += 1;
				}
			}
//TODO: BUG Exactly one is short and other - long!!!
			if ( counter == -diff ) {
				/*for ( auto edge = incomingEdges.begin(); edge != incomingEdges.end(); ++edge ) {
					if ( gp_.g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
						components.push_back(*edge);	
					}
				}*/
			}
		}
		if ( diff > 0 ) {
			int counter = 0;
			for ( auto edge = outgoing_edges.begin(); edge != outgoing_edges.end(); ++edge ) {
				if ( gp_.g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
					counter += 1;
				}
			}
			if ( counter == -diff ) {
				/*for ( auto edge = incomingEdges.begin(); edge != incomingEdges.end(); ++edge ) {
					if ( gp_.g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
						components.push_back(*edge);	
					}
				}*/
			}
		}
	}


	int CompareComponentLists( const vector<EdgeId>& v1, const vector<EdgeId>& v2 ) const {
		int counter = 0;
		for (auto it = v1.begin(); it != v1.end(); ++it) {
			if (find(v2.begin(), v2.end(), *it) == v2.end() && (quality_labeler_.quality(*it) > 0.5) ) {
				counter += 1;
				DEBUG(gp_.g.int_id(*it));
			}
		}
		DEBUG("\n");
		return counter;
	}

	public :
	template <class DetailedCoverage>
	void resolve_repeats_by_coverage( DetailedCoverage& coverage,
					size_t insert_size,
//TODO: Do not pass reference on graph_pack members.
//TODO: EdgeQuality to constructor
//TODO: check that everything works without reference:)

					EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
					PairedInfoIndexT<typename graph_pack::graph_t> & clustered_index,
					vector< PathInfo<typename graph_pack::graph_t> >& filtered_paths,
					string output_file_name) {

		auto filter = LoopFilter<graph_pack, DetailedCoverage>(gp_, coverage, tandem_lower_threshold_, tandem_upper_threshold_, repeat_length_upper_threshold_);
		filter.get_loopy_components( quality_labeler_ );
		INFO("Resolving repeats by coverage...");
		vector<EdgeId> components, singles;
		vector<EdgeId> components_ref, singles_ref;
		INFO("Getting components...");
		GetComponents(components, singles, labels_after, clustered_index );
//TODO:	CheckComponentsWithQuality()
//TODO:: Dima stopped here
		GetComponentsWithReference( components_ref, singles_ref);
		int fp(0), fn(0);
		DEBUG( "in components (size: " << components.size() << ") but not in components_ref:\n");
		fp = CompareComponentLists(components, components_ref);
		DEBUG("in components_ref (size: " << components_ref.size() << ") but not in components:\n");
		fn = CompareComponentLists(components_ref, components);
		DEBUG("False positives: " << (double) fp / components.size() << "False negatives: " << (double) fn / (fn + singles.size()) << "\n");
		//path with conjugate edges
		vector< vector<EdgeId> > resolved_paths;
		INFO("Processing graph...");
		TraverseComponents( components, singles, coverage, insert_size, resolved_paths );
		set<EdgeId> used_edges;
		INFO("Paths before joining:\n");
		for ( auto p = resolved_paths.begin(); p != resolved_paths.end(); ++p) {
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				INFO(gp_.g.int_id(*iter) << " (" << gp_.g.int_id(gp_.g.EdgeStart(*iter) ) << "," << gp_.g.int_id(gp_.g.EdgeEnd(*iter) ) << ") ");
			}
			INFO("\n");
		}
		INFO("Loops before joining:\n");
		for ( auto p = filter.resolved_loops().begin(); p != filter.resolved_loops().end(); ++p) {
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
			INFO(gp_.g.int_id(*iter) << "( " << gp_.g.length(*iter) << ") ");
			}
			INFO("\n");
		}
		vector< vector<EdgeId> > all_paths;
		JoinPaths(resolved_paths, filter.resolved_loops(), all_paths);
		FilterConjugate( used_edges, all_paths, filtered_paths);
		for ( auto p = all_paths.begin(); p != all_paths.end(); ++p) {
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				DEBUG(gp_.g.int_id(*iter) << " ");
			}
			DEBUG("\n");
		}
/*		INFO("-------------------------\n");
		for ( auto p = filtered_paths.begin(); p != filtered_paths.end(); ++p) {
			for ( auto iter = p->getPath().begin(); iter != p->getPath().end(); ++iter ) {
				INFO(gp_.g.int_id(*iter) << "( " << gp_.g.length(*iter) << "; " << coverage.GetInCov(*iter) << " " << coverage.GetOutCov(*iter) << ") ");
			}
			INFO("\n");
		}
	*/
		path_extend::PathContainer paths_output;
		for ( auto p = filtered_paths.begin(); p != filtered_paths.end(); ++p) {
			path_extend::BidirectionalPath* bidirectional_path = new path_extend::BidirectionalPath( gp_.g );
			path_extend::BidirectionalPath* conjugate_path = new path_extend::BidirectionalPath( gp_.g );
			auto tmpPath = p->getPath();
			for (auto it = tmpPath.begin(); it != tmpPath.end(); ++it ){
					bidirectional_path->PushBack(*it);
					EdgeId cedge = gp_.g.conjugate(*it);
					conjugate_path->PushFront( cedge );
			} 
			paths_output.AddPair( bidirectional_path, conjugate_path );
		}
		WriteResolved( used_edges, paths_output, output_file_name);
	}

	private:

	void filterConjugate( std::set<EdgeId>& usedEdges,
				const std::vector< std::vector<EdgeId> > & paths,
				std::vector< PathInfo<typename GraphPack::graph_t> >& filteredPaths) {


		GenomeConsistenceChecker<typename GraphPack::graph_t> checker(gp, 10, 0.2);
		INFO("filtering conjugate edges");
		for ( auto path = paths.begin(); path != paths.end(); ++path) {
			bool ifInsert = true;
			for (auto e = path->begin(); e != path->end(); ++e) {
				if ( gp_.g.conjugate(*e) == *e ) continue; 
				if ( used_edges.find(gp_.g.conjugate(*e)) != used_edges.end() ) {
//TODO:: this is not true, if we have autoreverse edge.					
					ifInsert = false;
					break;
				}
			}
			if (ifInsert) {
				INFO("inserting");
				if (!checker.IsConsistentWithGenome(*path)) {
					std::cout << "not consistent with genome: ";
					for (auto iter = path->begin(); iter != path->end(); ++iter) {
						auto positions = gp->edge_pos.GetEdgePositions(*iter);
						std::cout << gp->g.int_id(*iter) << " (";
						for (auto pos = positions.begin(); pos != positions.end(); ++pos) {
							std::cout << pos->mr.initial_range.start_pos << " - " << pos->mr.initial_range.end_pos << " ";
						}
						std::cout << ") ";
					}
					std::cout << std::endl;
				}

				filteredPaths.push_back(PathInfo<typename GraphPack::graph_t>(*path));
				for (auto e = path->begin(); e != path->end(); ++e)
					usedEdges.insert(*e);
				}
		}
		for ( auto path = filtered_paths.begin(); path != filtered_paths.end(); ++path ) {
			auto p = path->getPath();
			for (auto e = p.begin(); e != p.end(); ++e) {
				used_edges.insert(*e);
			}
		}
	}

	void GetComponentsWithReference( vector<EdgeId>& components, vector<EdgeId>& singles){
		DEBUG("Finding Components With Paired Info");
		for (auto iter = gp_.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (quality_labeler_.quality(*iter) > 1.5 ) {
				components.push_back(*iter);
			}
			else {
				singles.push_back(*iter);
			}
		}
	}

	bool IsBulge( const EdgeId& edge ){
		auto edgeStart = gp_.g.EdgeStart(edge);
		auto edgeEnd = gp_.g.EdgeEnd(edge);
		auto outgoingFromStart = gp_.g.OutgoingEdges(edgeStart);
		for (auto e = outgoingFromStart.begin(); e != outgoingFromStart.end(); ++e){
			if (*e == edge) continue;
			if (gp_.g.EdgeEnd(*e) == edgeEnd) return true;

		}

		return false;
	}

	template< class Graph>
	bool checkIfComponentByPairedInfo( EdgeId edge, PairedInfoIndexT<Graph>& clustered_index, std::set<EdgeId>& prohibitedEdges ) {

//		auto improver = PairInfoImprover<Graph>(gp->g, clustered_index);
		auto inner_map = clustered_index.GetEdgeInfo(edge, 0);
		for (auto I_1 = inner_map.Begin(), E = inner_map.End(); I_1 != E; ++I_1) {
			for (auto I_2 = inner_map.Begin(); I_2 != E; ++I_2) {
				if (I_1 == I_2) continue;
				EdgeId e1 = (*I_1).first;
				const Point& p1 = (*I_1).second;
				EdgeId e2 = (*I_2).first;
				const Point& p2 = (*I_2).second;

				if (prohibitedEdges.find(e1) != prohibitedEdges.end() || prohibitedEdges.find(e2) != prohibitedEdges.end() ) continue;
				if ( p1.d * p2.d < 0 || p2.d > p1.d ) continue;
//				if (!improver.IsConsistent(edge, e1, e2, p1, p2)) {
//					std::cout << "Inconsistent for " << gp->g.int_id(edge) << ": " << gp->g.int_id(e1) << " " << gp->g.int_id(e2) << std::endl;
//					return true;
//				}
			}
		}
		return false;
	}

	void GetComponents( vector<EdgeId>& components, vector<EdgeId>& singles, 
				EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
				PairedInfoIndexT<Graph>& clustered_index){

		typedef int times;
		map<VertexId, times> out_degree;
		map<VertexId, times> in_degree;
		for (auto e_iter = gp_.g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
			VertexId from = gp_.g.EdgeStart(*e_iter);
			VertexId into = gp_.g.EdgeEnd(*e_iter);
			if (out_degree.find(from) != out_degree.end()){
				out_degree[from] += 1;
			}
			else{
				out_degree[from] = 1;
			}
			if (in_degree.find(into) != in_degree.end()){
				in_degree[into] += 1;
			}
			else{
				in_degree[into] = 1;
			}

		}
		int byPairedInfo = 0;
		int byTopology = 0;
		int confirmedByPairedInfo = 0;
		for (auto e_iter = gp_.g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
			VertexId from = gp_.g.EdgeStart(*e_iter);
			VertexId into = gp_.g.EdgeEnd(*e_iter);
			if ( gp_.g.length(*e_iter) >= cfg::get().rr.max_repeat_length ) {
				singles.push_back(*e_iter);
			}

			else if (! ( ((in_degree.find(from) == in_degree.end()) || out_degree[from] > 1) && ((out_degree.find(into) == out_degree.end()) || in_degree[into] > 1) ) ) {
				byTopology += 1;
				edge_to_kind_[*e_iter] = TOPOLOGY;
				components.push_back(*e_iter);
			}
			else if ( IfRepeatByPairedInfo(*e_iter, clustered_index ) ) {
				components.push_back(*e_iter);
				byPairedInfo += 1;
				edge_to_kind_[*e_iter] = PAIREDINFO;
			}
			else if( gp_.g.length(*e_iter) < repeat_length_upper_threshold_ && (in_degree.find(from) != in_degree.end()) && (out_degree.find(into) != out_degree.end()) ){
				edge_to_kind_[*e_iter] = LENGTH;
				components.push_back(*e_iter);
			} 
			else{
				singles.push_back(*e_iter);
			}

		}
		DEBUG("Number of edges identified by paired info: " << byPairedInfo);
		DEBUG("Number of edges identified by topology: " << byTopology);
		DEBUG("Number of edges identified by topology confirmed by paired info: " << confirmedByPairedInfo);
		DEBUG("SINGLES: ");
		for (auto sit = singles.begin(); sit != singles.end(); ++sit){
			DEBUG(gp_.g.int_id(*sit) << ", ");
		}
		DEBUG("COMPONENTS: ");
		for (auto cit = components.begin(); cit != components.end(); ++cit){
			DEBUG(gp_.g.int_id(*cit) << ", ");
		}
		DEBUG("\n");
	}

	template <typename T1, typename T2>
	struct CompareSecond {
		typedef pair<T1, T2> type;
		bool operator ()(type const& a, type const& b) const {
			return a.second < b.second;
	    	}
	};

	// idea suggested by Anton K: count the probablities of distribution of all the possible ways of repeat resolution for a given repeat 		
	// look if the function is bended left (goof for repeat resolution) or right (bad)
	void CountProbabililtiesDistribution ( FILE* file, const vector< vector <double> > transition_probabilities ) {
		unsigned size = transition_probabilities.size();
		vector<unsigned> permutation;
		for ( unsigned i = 0; i != size; ++i ) {
			permutation.push_back(i);
		}
		vector<double> probabilities;   
		do {
			probabilities.push_back(1);
			for ( unsigned i = 0; i != size; ++i ) {
				probabilities[probabilities.size()-1] *= transition_probabilities[permutation[i]][i] ;	
			}
			
		}
		while ( std::next_permutation(permutation.begin(), permutation.end() ));
		DEBUG("\n");
		sort(probabilities.begin(), probabilities.end());

			VertexId from = gp->g.EdgeStart(*e_iter);
			VertexId into = gp->g.EdgeEnd(*e_iter);

			//std::cout << e_iter->int_id() << std::endl;
			if ( gp->g.length(*e_iter) >= cfg::get().max_repeat_length ) {
				singles.push_back(*e_iter);
			}
			pred = *it;
		}
		DEBUG("\n");
	}

/*
	void ChooseMostLikelyPairs( const vector< vector <double> > transition_probabilities,
					 vector<pair<EdgeId,EdgeId>>& pairs_of_edges, 
					 const vector<EdgeId>& incoming_edges,
					 const vector<EdgeId>& outgoing_edges) 
	{

		unsigned k = 0;
		set<unsigned> matched;
		cout << "MATCH: " << endl;
		for ( auto vec = transition_probabilities.begin(); vec != transition_probabilities.end(); ++vec, ++k ) {
		
			double max_prob = (*vec)[0];
			unsigned max_index = 0;
			for (unsigned i = 1; i < vec->size(); ++i) {
				if (max_prob < (*vec)[i]) {
					max_prob = (*vec)[i];
					max_index = i;
				}
			}

			bool if_matches = true;

			if (max_prob < 0.15) if_matches = false;

			// actually must check if the new value of max prob is greater than in previous match
			if (matched.find(max_index) != matched.end()) if_matches = false;
			if (if_matches)
			for (unsigned i = 0; i < vec->size(); ++i) {
				if (i == max_index) continue;
				if ( abs(max_prob - (*vec)[i]) < 0.1 )
					if_matches = false;
			}

			if (if_matches)
			for (unsigned j = 0; j < vec->size(); ++j) {
				if (j == k) continue;
				if ( abs(max_prob - transition_probabilities[j][max_index]) < 0.1 )
					if_matches = false;
			}

			if (if_matches) {
				cout << k << " " <<  max_index << endl;
				matched.insert(k);
				matched.insert(max_index);
				pairs_of_edges.push_back(make_pair(incoming_edges[k], outgoing_edges[max_index]) );
			}

		}
	}



	template <class DetailedCoverage>
	bool MatchPairs( vector<EdgeId>& incoming_edges,
			 vector<EdgeId>& outgoing_edges,
			 vector<pair<EdgeId,EdgeId>>& pairs_of_edges,
			 const vector<EdgeId>& component,
			 BucketMapper<Graph> &bm,
			 DetailedCoverage& coverage,
			 EdgeQuality<typename graph_pack::graph_t>& quality_labeler,
			 FILE* file)  {


//		if (incomingEdges.size() > 5 || outgoingEdges.size() > 5) return false;
		double shift = 25;

		vector< vector <double> > transition_probabilities ;
		for ( unsigned i = 0; i < incoming_edges.size(); i++) {
			 transition_probabilities.push_back(vector<double>(outgoing_edges.size(),1));
		}

		int in_edge_counter(0);


		//cout << "component size: " << component.size() << endl;
		for ( auto in_edge = incoming_edges.begin(); in_edge != incoming_edges.end(); ++in_edge, ++in_edge_counter ) {
	
			//cout << "incoming edge " << gp_.g.int_id(*in_edge) << endl;
	
			double in_cov = coverage.GetOutCov(*in_edge);
			int in_bucket = bm.GetCoverageBucket(in_cov);
			 
			 	int out_edge_counter(0);
				for ( auto out_edge = outgoing_edges.begin(); out_edge != outgoing_edges.end(); ++out_edge, ++out_edge_counter ) {
				
						
						//cout << "outgoing edge " << gp_.g.int_id(*out_edge) << endl;
						double out_cov = coverage.GetInCov(*out_edge);
						//cout << "out_cov: " << out_cov << endl;
						int out_bucket = bm.GetCoverageBucket(out_cov);
						//cout << "out_bucket: " << out_bucket << endl;

						int distance(0);
						CountDistance(*in_edge, *out_edge, component, distance);

						//cout << distance << " " << in_bucket << " " << out_bucket << endl;
						double probability = bm.GetProbabilityFromBucketToBucketForDistance (in_bucket, out_bucket, distance, shift) ;
						//cout << probability << endl;
						transition_probabilities[in_edge_counter][out_edge_counter] = probability;
		        	} 
				//cout << endl;
		}

		CountProbabililtiesDistribution ( file, transition_probabilities );
		GetComponentInfo(component, incoming_edges, outgoing_edges, coverage, transition_probabilities, quality_labeler );

		ChooseMostLikelyPairs( transition_probabilities, pairs_of_edges, incoming_edges, outgoing_edges);

		if (pairs_of_edges.size() > 0) return true;

		return false;
		
	}


	void findClosest(vector<pair<EdgeId, coverage_value>>& incomingEdgesCoverage,
			vector<pair<EdgeId, coverage_value>>& outgoingEdgesCoverage,
			vector<pair<EdgeId,EdgeId>>& pairsOfEdges){

		int Length = (int) min(incomingEdgesCoverage.size(), outgoingEdgesCoverage.size());

		//TODO: Move to config
		//double threshold_one_list_(0.80), threshold_match_(0.64);
		//double threshold_one_list_(0.80), threshold_match_(0.70);

		for (int i = 0; i < Length - 1; ++i) {

			double valueOneListIn = min (incomingEdgesCoverage[i].second, incomingEdgesCoverage[i+1].second) / max (incomingEdgesCoverage[i].second, incomingEdgesCoverage[i+1].second) ;
			double valueOneListOut = min (outgoingEdgesCoverage[i].second, outgoingEdgesCoverage[i+1].second) / max (outgoingEdgesCoverage[i].second, outgoingEdgesCoverage[i+1].second) ;
			double valuePair = min (incomingEdgesCoverage[i].second, outgoingEdgesCoverage[i].second) / max (incomingEdgesCoverage[i].second, outgoingEdgesCoverage[i].second) ;

			if ( valueOneListIn > threshold_one_list_ || valueOneListOut > threshold_one_list_ || valuePair < threshold_match_ ){
				return;
			}
		}
		int i = Length - 1;
		double valuePair = min (incomingEdgesCoverage[i].second, outgoingEdgesCoverage[i].second) / max (incomingEdgesCoverage[i].second, outgoingEdgesCoverage[i].second) ;
		if ( valuePair < threshold_match_ ){
			return;
		}

		if ( Length == (int) incomingEdgesCoverage.size() ){
			for (unsigned i = Length - 1; i < outgoingEdgesCoverage.size() - 1; ++i){
				double valueOneList = (double) min (outgoingEdgesCoverage[i].second, outgoingEdgesCoverage[i+1].second) / max (outgoingEdgesCoverage[i].second, outgoingEdgesCoverage[i+1].second);
				if (valueOneList > threshold_one_list_){
					return;
				}
			}
		}
		else {

			for (unsigned i = Length - 1; i < incomingEdgesCoverage.size() - 1; ++i){
				double valueOneList = (double) min (incomingEdgesCoverage[i].second, incomingEdgesCoverage[i+1].second) / max (incomingEdgesCoverage[i].second, incomingEdgesCoverage[i+1].second);
				if (valueOneList > threshold_one_list_){
					return;
				}
			}


		}

		for (int i = 0; i < Length; ++i){
			pairsOfEdges.push_back( make_pair(incomingEdgesCoverage[i].first, outgoingEdgesCoverage[i].first)  );
		}

	}


	bool ContainsSmallLoop( const vector<EdgeId>& path){
		
		if ( path.size() == 1 ) return true;

		if ( path.size() == 2 ) {

			if (gp_.g.EdgeStart(path[0]) == gp_.g.EdgeEnd(path[1]) && gp_.g.EdgeStart(path[1]) == gp_.g.EdgeEnd(path[0]))

				return true;
		}

		return false;
	}

	bool ContainsOnlyShortEdges( const vector<EdgeId>& path){

		for ( auto it = path.begin(); it != path.end(); ++it ) {
		    if (gp_.g.length(*it) >= repeat_length_upper_threshold_ || edge_to_kind_[*it] == TOPOLOGY )
		        return false;
			
		}
		return true;
	}


	void bfs ( const EdgeId& edge,  set<EdgeId>& visited_edges, const vector<EdgeId>& component, int& curLen, int& maxPathLen) {

		visited_edges.insert(edge);
		auto incoming_edges = gp_.g.IncomingEdges(gp_.g.EdgeStart(edge));

		for ( auto e = incoming_edges.begin(); e != incoming_edges.end(); ++e) {
		    if ( find(component.begin(), component.end(), *e) != component.end() && visited_edges.find(*e) == visited_edges.end() ){
				curLen += gp_.g.length(*e);
				if (curLen > maxPathLen) maxPathLen = curLen;
				bfs(*e, visited_edges, component, curLen, maxPathLen);
			}

		}

		auto outgoingEdges = gp_.g.OutgoingEdges(gp_.g.EdgeEnd(edge));
		for ( auto e = outgoingEdges.begin(); e != outgoingEdges.end(); ++e) {
	
			if ( find(component.begin(), component.end(), *e) != component.end() && visited_edges.find(*e) == visited_edges.end() ){
					
					curLen += gp_.g.length(*e);
					if (curLen > maxPathLen) maxPathLen = curLen;
					bfs(*e, visited_edges, component, curLen, maxPathLen);
			}
		}

	}

	int GetLongestPathLength( const vector<EdgeId>& component ){
	// gets a repetitive component and calculates the length of the longest path in it

		set<EdgeId> visited_edges;
		vector<vector<EdgeId>> paths;
		int maxPathLen = gp_.g.length(component[0]);
		for ( auto edge = component.begin(); edge != component.end(); ++edge ){

			if (visited_edges.find(*edge) != visited_edges.end()) continue;
			int curLen = gp_.g.length(*edge);
			visited_edges.insert(*edge);
			bfs(*edge, visited_edges, component, curLen, maxPathLen);

		}

		return maxPathLen;

	}


	bool CheckRepeatDetection( const vector<EdgeId>& component, 
				const vector<EdgeId>& incoming_edges, 
				const vector<EdgeId>& outgoing_edges,
				const EdgeQuality<typename graph_pack::graph_t>& quality_labeler ) {

		for ( auto e = component.begin(); e != component.end(); ++e ) {

			if ( quality_labeler.quality(*e) <= 1.5 ) {
				return false;
			}
		
		}

		for ( auto e = incoming_edges.begin(); e != incoming_edges.end(); ++e ) {
			
			if ( quality_labeler.quality(*e) > 1.5 ) {
				return false;
			}

		}

		for ( auto e = outgoing_edges.begin(); e != outgoing_edges.end(); ++e ) {
			
			if ( quality_labeler.quality(*e) > 1.5 ) {
				return false;
			}
		}
		return true;
	}


	template <class DetailedCoverage>
	void GetComponentInfo(const vector<EdgeId>& component, const vector<EdgeId>& incoming_edges, const vector<EdgeId>& outgoing_edges, const DetailedCoverage& coverage, 
		const vector< vector <double> >& transition_probabilities, const EdgeQuality<typename graph_pack::graph_t>& quality_labeler ) {

			cout << "Component: " << endl;
			for ( auto iter = component.begin(); iter != component.end(); ++iter ) {
				cout << gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) << 
							" average edge coverage " << gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << " "; 
				
				if (edge_to_kind_[*iter] == TOPOLOGY) 
					cout << "TOPOLOGY" << endl;
				else if  (edge_to_kind_[*iter] == LENGTH )
					cout << "LENGTH" << endl;
				else if (edge_to_kind_[*iter] == PAIREDINFO )
					cout << "PAIREDINFO" << endl;
		}
			cout << endl;
			cout << "incoming edges: " << endl;
			for ( auto iter = incoming_edges.begin(); iter != incoming_edges.end(); ++iter ) {
				cout << gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) << " outgoing edge coverage: " << coverage.GetOutCov(*iter) << 
							" average edge coverage " << gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << endl;
			}
			cout << endl;
			cout << "outgoing edges: " << endl;
			for ( auto iter = outgoing_edges.begin(); iter != outgoing_edges.end(); ++iter ) {
				cout << gp_.g.int_id(*iter)  << " edge length: " << gp_.g.length(*iter) << " incoming edge coverage: " << coverage.GetInCov(*iter) << 
							" average edge coverage " << gp_.g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << endl;
			}
			cout << endl;
			
			bool correct_component = CheckRepeatDetection( component, incoming_edges, outgoing_edges, quality_labeler  );

			if (!correct_component) {
				cout << "repeat is detected incorrectly" << endl;
			}
			if (transition_probabilities.size() > 0) {
				cout << "transition probabilities" << incoming_edges.size() << "x" << outgoing_edges.size() << ":" << endl;
				for (auto vec = transition_probabilities.begin(); vec != transition_probabilities.end(); ++vec ) {
			
					for ( auto prob = vec->begin(); prob != vec->end(); ++prob ) {

						printf("%4.5f ", *prob);
					}
					cout << endl;
				}
			}
			cout << endl;

	}

*/

	//todo WTF?!!!
	template <class DetailedCoverage>
	void TraverseComponents( const vector<EdgeId>& components, 
				const vector<EdgeId>& singles, const DetailedCoverage& coverage, 
				size_t insert_size,
				vector< vector<EdgeId>> & resolved_paths ) {

		set<EdgeId> visited_edges;
		DEBUG("Traversing components");
		FILE *file = fopen((cfg::get().output_dir+"repeats.log").c_str(), "w");
		//FILE* file = fopen("/home/ksenia/path_resolved.log", "w");
		//FILE* file = fopen("/home/ksenia/probabilities_22.log", "w");
        	int number_of_buckets = 20;
		int K_ = (int) cfg::get().K + 1;
		BucketMapper<conj_graph_pack::graph_t, KmerIndex> bm (gp_.g, kmer_index_, K_, number_of_buckets);
		bm.InitBuckets();
		int pure_tandem(0), repetitive_tandem(0), ordinal_repeat(0);
		for ( auto edge = components.begin(); edge != components.end(); ++edge ) {
			if ( visited_edges.find(*edge) != visited_edges.end() ){
				continue;
			}
			vector<EdgeId> incoming_edges, outgoing_edges;
			vector<EdgeId> path;
			set<VertexId> component_vertices;
			bool if_loop = false;
			visit(*edge, visited_edges, component_vertices, path, components, singles, incoming_edges, outgoing_edges, if_loop);
			Repeat<graph_pack> repeat(gp_, incoming_edges, outgoing_edges, path, repeat_length_upper_threshold_, edge_to_kind_,file);
			vector<vector <double>> transition_probabilities ;
			if (if_loop) {
				if (incoming_edges.size() == (size_t) 1 && incoming_edges.size() == outgoing_edges.size() )
					pure_tandem += 1;
				else repetitive_tandem += 1;
				DEBUG("loop!\n");
				continue;
			}
			if (path.size() == (size_t) 0 ) continue;
			if ( incoming_edges.size() < (size_t) 2 || outgoing_edges.size() < (int) 2) continue;
			ordinal_repeat += 1;
			if ( repeat.IfContainsOnlyShortEdges() ) {
				DEBUG("contains only short edges");
				DEBUG("component: ");
				for ( auto iter = path.begin(); iter != path.end(); ++iter ) {
					DEBUG(gp_.g.int_id(*iter)  << " ");
				}
				DEBUG("\n");
				continue;
			}
			if (incoming_edges.size() != outgoing_edges.size()){
				continue;
			}
			//vector<vector<EdgeId>> resolved_paths;
			repeat.Resolve(bm, coverage, quality_labeler_, resolved_paths);
		}
		
		DEBUG( "pure tandems: " << pure_tandem << "\nrepeats + tandems: " << repetitive_tandem << "\nordinal repeats: " << ordinal_repeat << "\n");
		fclose(file);
	}

	//private:
	//DECL_LOGGER("CoverageBasedRR");
/*
	// verify if the sequence of edges is genomic path
	bool MatchReference( const path_extend::BidirectionalPath& path) {

		auto ref_pos = gp_->edge_pos;
		EdgeId edge = path.At(0);
		auto pos_it = ref_pos.edges_positions().find(edge);
		//VERIFY(pos_it != ref_pos.edges_positions().end();

		if ( pos_it->second.size() == 1 ){

			auto nextStart = pos_it->second[0].end();
			return match( path, 1, nextStart + 1, ref_pos );
		}
		//INFO("first is repeat - fail or uncovered");
		cout << pos_it->second.size() << endl;
		for (size_t i = 0; i < pos_it->second.size(); i++) {
		      cout << "    " << pos_it->second[i].contigId_ << " "
		    			<< pos_it->second[i].start() << " - "
				         << pos_it->second[i].end() << endl;
		}
		return false;
	}
*/
	};
}
#endif
