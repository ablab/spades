#ifndef COVERAGE_BASED_RR
#define COVERAGE_BASED_RR

#include <vector>
#include <map>
#include <set>
#include <list>
#include "graph_print_utils.hpp"
#include "indices/debruijn_kmer_index.hpp"
#include "graph_pack.hpp"
#include "pair_info_improver.hpp"
#include "path_extend/pe_io.hpp"
#include "path_extend/bidirectional_path.hpp"
#include "graphio.hpp"
#include "long_read_storage.hpp"
#include "bucket_mapper.hpp"

namespace debruijn_graph{

template <class graph_pack>
class CoverageBasedResolution {
	typedef double coverage_value;
	typedef enum { TOPOLOGY, LENGTH, PAIREDINFO } kind_of_repeat; 

	//TODO: remove 
	const graph_pack& gp_;
	
	const DeBruijnEdgeIndex<EdgeId>& kmer_index_;

	map<EdgeId, kind_of_repeat> edge_to_kind_;

	const double threshold_one_list_;
	const double threshold_match_;
	const double threshold_global_;
	const double tandem_lower_threshold_;
	const double tandem_upper_threshold_;
	const double repeat_length_upper_threshold_;
	
	public:
	CoverageBasedResolution( const graph_pack& gpack_arg, const DeBruijnEdgeIndex<EdgeId>& kmer_index, double threshold_one_list, double threshold_match, 
				double threshold_global, double tandem_lower_threshold, double tandem_upper_threshold, double repeat_length_upper_threshold) :

											gp_(gpack_arg),
											kmer_index_(kmer_index),
											threshold_one_list_(threshold_one_list), 
											threshold_match_(threshold_match),
											threshold_global_(threshold_global), 
											tandem_lower_threshold_(tandem_lower_threshold), 
											tandem_upper_threshold_(tandem_upper_threshold),
											repeat_length_upper_threshold_(repeat_length_upper_threshold){
	}
	
	
	private:
	void JoinPaths( vector<vector<EdgeId> >& paths, //vector< vector<EdgeId> >& paths,
			vector< vector<EdgeId> >& resolvedLoops,
			vector< vector<EdgeId> >& all_paths ) {

		map< EdgeId, vector<EdgeId> > startEdgeToPath;
		map< EdgeId, vector<EdgeId> > backEdgeToPath;

		vector< vector<EdgeId> > bothPaths;
		for (auto path = resolvedLoops.begin(); path != resolvedLoops.end(); ++path) {
			bothPaths.push_back(*path);

		}
		for (auto path = paths.begin(); path != paths.end(); ++path) {
			bothPaths.push_back(*path);
		}

		INFO("before map filling");


		for (auto path = bothPaths.begin(); path != bothPaths.end(); ++path) {
		//	cout << gp_.g.int_id(path->front()) << endl;
			//VERIFY(startEdgeToPath.find(path->front()) == startEdgeToPath.end());
			startEdgeToPath[path->front()] = *path;
		}

		for (auto path = bothPaths.begin(); path != bothPaths.end(); ++path) {
			//VERIFY(backEdgeToPath.find(path->back()) == backEdgeToPath.end());
			backEdgeToPath[path->back()] = *path;
		}

		for (auto path = bothPaths.begin(); path != bothPaths.end(); ++path) {
			if ( backEdgeToPath.find( path->front() ) != backEdgeToPath.end()) {
				continue;
			}

			/*INFO("path before");
			for ( auto e = path->begin(); e != path->end(); ++e ){
				cout << gp_.g.int_id(*e) <<  " ";
			}
			cout << endl; */
			bool updated = true;
			while (updated) {
				auto foundPath = startEdgeToPath.find( path->back() );
				if (foundPath != startEdgeToPath.end() ) {
					//INFO("found path before");
					/*for ( auto e = foundPath->second.begin(); e != foundPath->second.end(); ++e ){
						cout << gp_.g.int_id(*e) <<  " ";
					}
					cout << endl;*/
					path->insert(path->end(), foundPath->second.begin() + 1, foundPath->second.end());
				} else {
					updated = false;
				}
			}
			/*INFO("path after");
			for ( auto e = path->begin(); e != path->end(); ++e ){
				cout << gp_.g.int_id(*e) <<  " ";
			}
			cout << endl;*/
			all_paths.push_back(*path);
		}
		INFO("out of path joining");
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
				cout << "INCOMING COMPONENT UPDATED" << endl;
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
				cout << "OUTGOING COMPONENT UPDATED" << endl;
				/*for ( auto edge = incomingEdges.begin(); edge != incomingEdges.end(); ++edge ) {

					if ( gp_.g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
						components.push_back(*edge);	
					}
				}*/
			}
		}
	}


	public :
	template <class DetailedCoverage>
	void resolve_repeats_by_coverage( DetailedCoverage& coverage, 
					size_t insert_size,
//TODO: Do not pass reference on graph_pack members.
//TODO: EdgeQuality to constructor
//TODO: check that everything works without reference:)

					EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
					EdgeQuality<typename graph_pack::graph_t>& quality_labeler,
					PairedInfoIndexT<typename graph_pack::graph_t> & clustered_index,

					vector< PathInfo<typename graph_pack::graph_t> >& filtered_paths
					//set<EdgeId>& prohibitedEdges,
					//vector< vector<EdgeId>>& resolvedLoops ) 
					) {

		auto filter = LoopFilter<graph_pack, DetailedCoverage>(gp_, coverage, tandem_lower_threshold_, tandem_upper_threshold_, repeat_length_upper_threshold_);
		filter.get_loopy_components( quality_labeler );

		INFO("Resolving repeats by coverage...");
		vector<EdgeId> components, singles;
		vector<EdgeId> componentsRef, singlesRef;
		INFO("Getting components...");

		GetComponents(components, singles, labels_after, quality_labeler, clustered_index, filter.prohibitedEdges );
//TODO:	CheckComponentsWithQuality()
//TODO:: Dima stopped here
		GetComponentsWithReference( componentsRef, singlesRef, quality_labeler, filter.prohibitedEdges );

		int fp(0), fn(0);

		cout << "in components (size: " << components.size() << ") but not in componentsRef: " << endl;
		for (auto it = components.begin(); it != components.end(); ++it) {

			if (find(componentsRef.begin(), componentsRef.end(), *it) == componentsRef.end() 
				&& (quality_labeler.quality(*it) > 0.5) ) {
				fp += 1;
				cout << gp_.g.int_id(*it);
				if ( edge_to_kind_[*it] == TOPOLOGY ) {
					  cout << " (TOPOLOGY) , ";
				}
				if ( edge_to_kind_[*it] == LENGTH ) {
					cout << " (LENGTH) , ";
				}

			}

		}
		cout << endl;
		
	
		cout << "in componentsRef (size: " << componentsRef.size() << ") but not in components: " << endl;
		for (auto it = componentsRef.begin(); it != componentsRef.end(); ++it) {

			if (find(components.begin(), components.end(), *it) == components.end()) {
				fn += 1;
				cout << gp_.g.int_id(*it) << ", ";
			}

		}
		cout << endl;
		cout << "False positives: " << (double) fp / components.size() << "False negatives: " << (double) fn / (fn + singles.size()) << endl;




		//path with conjugate edges
		vector< vector<EdgeId> > resolved_paths;
	
		//getComponents( gp_, components, singles, quality_labeler, unresolvedLoops );

		INFO("Traversing graph...");
		TraverseComponents( components, singles, coverage, insert_size, resolved_paths, quality_labeler );

		set<EdgeId> used_edges;

		cout << "Paths before joining: " << endl;
		for ( auto p = resolved_paths.begin(); p != resolved_paths.end(); ++p) {
			//fprintf(file, "resolved path \n");
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				cout << gp_.g.int_id(*iter) << " (" << gp_.g.int_id(gp_.g.EdgeStart(*iter) ) << "," << gp_.g.int_id(gp_.g.EdgeEnd(*iter) ) << ") ";
				//fprintf(file, "%d ", gp_.g.int_id(*iter));
				//fprintf(file, " ");
			}
			cout << endl;
			//fprintf(file,"\n");
		}
		cout << "Loops before joining: " << endl;
		for ( auto p = filter.resolvedLoops.begin(); p != filter.resolvedLoops.end(); ++p) {
			//fprintf(file, "resolved path \n");
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				cout << gp_.g.int_id(*iter) << "( " << gp_.g.length(*iter) << ") ";
				//fprintf(file, "%d ", gp_.g.int_id(*iter));
				//fprintf(file, " ");
			}
			cout << endl;
			//fprintf(file,"\n");
		}

		vector< vector<EdgeId> > all_paths;
		JoinPaths(resolved_paths, filter.resolvedLoops, all_paths);
		FilterConjugate( used_edges, all_paths, filtered_paths);
		//cout << "before filtering size " << allPaths.size() << " filtered size: " << filteredPaths.size() << endl;

		
		for ( auto p = all_paths.begin(); p != all_paths.end(); ++p) {
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				cout << gp_.g.int_id(*iter) << " ";
			}
			cout << endl;
		}

		cout << "-------------------------" << endl;
		for ( auto p = filtered_paths.begin(); p != filtered_paths.end(); ++p) {
			//fprintf(file, "resolved path \n");
			for ( auto iter = p->getPath().begin(); iter != p->getPath().end(); ++iter ) {
				cout << gp_.g.int_id(*iter) << "( " << gp_.g.length(*iter) << "; " << coverage.GetInCov(*iter) << " " << coverage.GetOutCov(*iter) << ") ";
				//fprintf(file, "%d ", gp_.g.int_id(*iter));
				//fprintf(file, " ");
			}
			cout << endl;
			//fprintf(file,"\n");
		}

		//fclose(file);
		string file_name = cfg::get().output_dir + "resolved_by_coverage.fasta";
		//INFO("Writing result");

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

		WriteResolved( used_edges, paths_output, file_name);
	}

	private:

	void FilterConjugate( set<EdgeId>& used_edges,
				const vector< vector<EdgeId> > & paths,
				vector< PathInfo<typename graph_pack::graph_t> >& filtered_paths) {


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
				/*if (! gp_->edge_pos.IsConsistentWithGenome(*path)) {
					for (auto iter = path->begin(); iter != path->end(); ++iter) {
						auto positions = gp_->edge_pos.GetEdgePositions(*iter);
				}*/

				filtered_paths.push_back(PathInfo<typename graph_pack::graph_t>(*path));
				for (auto e = path->begin(); e != path->end(); ++e) 
					used_edges.insert(*e);
				}
		}
		
		//INFO("inserting paths into set of the used edges");
		for ( auto path = filtered_paths.begin(); path != filtered_paths.end(); ++path ) {
			auto p = path->getPath();
			//cout << "size of p: " << p.size() << endl;
			for (auto e = p.begin(); e != p.end(); ++e) {
				used_edges.insert(*e);
			}
		}

		//INFO("out of filtering");
	}

	void GetComponentsWithReference( vector<EdgeId>& components, vector<EdgeId>& singles,
					EdgeQuality<typename graph_pack::graph_t>& quality_labeler,
					set<EdgeId>& prohibitedEdges ){

		INFO("Finding Components With Paired Info");
		for (auto iter = gp_.g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {

			if ( prohibitedEdges.find(*iter) != prohibitedEdges.end() ){
				continue;
			}

			if (quality_labeler.quality(*iter) > 1.5 ) {

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
		//auto incomingToEnd = gp_.g.IncomingEdges(edgeEnd);
 
		for (auto e = outgoingFromStart.begin(); e != outgoingFromStart.end(); ++e){
			if (*e == edge) continue;
			if (gp_.g.EdgeEnd(*e) == edgeEnd) return true;

		}

		return false;
	}

	template< class Graph>
	bool CheckIfComponentByPairedInfo( EdgeId edge, PairedInfoIndexT<Graph>& clustered_index, set<EdgeId>& prohibitedEdges ) {

		io::SequencingLibrary<debruijn_config::DataSetData> lib;
		auto improver = PairInfoImprover<Graph>(gp_.g, clustered_index,lib);
		InnerMap<Graph> inner_map = clustered_index.GetEdgeInfo(edge, 0);
		for (auto I_1 = inner_map.Begin(), E = inner_map.End(); I_1 != E; ++I_1) {
			for (auto I_2 = inner_map.Begin(); I_2 != E; ++I_2) {
				if (I_1 == I_2) continue;
				EdgeId e1 = (*I_1).first;
				const Point& p1 = (*I_1).second;
				EdgeId e2 = (*I_2).first;
				const Point& p2 = (*I_2).second;
				
				if (prohibitedEdges.find(e1) != prohibitedEdges.end() || prohibitedEdges.find(e2) != prohibitedEdges.end() ) continue;
				if ( p1.d * p2.d < 0 || p2.d > p1.d ) continue;
				if (!improver.IsConsistent(edge, e1, e2, p1, p2)) {
					//cout << "Inconsistent for " << gp_.g.int_id(edge) << ": " << gp_.g.int_id(e1) << " " << gp_.g.int_id(e2) << endl;
					return true;
				}
			}
		}

		return false;
	}

	void GetComponents( vector<EdgeId>& components, vector<EdgeId>& singles, 
				EdgeLabelHandler<typename graph_pack::graph_t>& labels_after,
				EdgeQuality<Graph>& quality_labeler,
				PairedInfoIndexT<Graph>& clustered_index,
				set<EdgeId>& prohibitedEdges ){

		typedef int times;
		map<VertexId, times> out_degree;
		map<VertexId, times> in_degree;
	/*	double LengthCutoff = 0.0;

		int numEdges = 0;
		for (auto e_iter = gp_.g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
			LengthCutoff += gp_.g.length(*e_iter);	
			numEdges += 1;
		}

		LengthCutoff = 0.25 * LengthCutoff / numEdges;	
	*/

//		cout << "Length Cutoff: "  << LengthCutoff << endl;
		INFO("Getting Components");
		INFO("Counting degrees of vertices");
		for (auto e_iter = gp_.g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
		
			//if ( prohibitedEdges.find(*e_iter) != prohibitedEdges.end() ) continue;
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

		/*int numEdges = 0;
		for (auto e_iter = gp_.g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
			LengthCutoff += gp_.g.length(*e_iter);	
			numEdges += 1;
		}

		cout << "Length cutoff: " << 0.25 * LengthCutoff / numEdges << endl;
		*/
		int byPairedInfo = 0;
		int byTopology = 0;
		int confirmedByPairedInfo = 0;
		for (auto e_iter = gp_.g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){

			if ( prohibitedEdges.find(*e_iter) != prohibitedEdges.end() ){ 
				continue;
			}

			VertexId from = gp_.g.EdgeStart(*e_iter);
			VertexId into = gp_.g.EdgeEnd(*e_iter);

			if ( gp_.g.length(*e_iter) >= cfg::get().rr.max_repeat_length ) {

				singles.push_back(*e_iter);
			}
			else if (! ( ((in_degree.find(from) == in_degree.end()) || out_degree[from] > 1) && ((out_degree.find(into) == out_degree.end()) || in_degree[into] > 1) ) ) {
				byTopology += 1;
				//if (checkIfComponentByPairedInfo(*e_iter, clustered_index, prohibitedEdges )) confirmedByPairedInfo++; 
				edge_to_kind_[*e_iter] = TOPOLOGY;
				components.push_back(*e_iter);
			}
			/*else if ( checkIfComponentByPairedInfo(*e_iter, clustered_index, prohibitedEdges ) ) {
				components.push_back(*e_iter);
				//	cout << "Component Edge Detected By Paired Info: " << gp_.g.int_id(*e_iter) << endl;
				byPairedInfo += 1;
				edge_to_kind_[*e_iter] = PAIREDINFO;
			}*/

			else if( gp_.g.length(*e_iter) < /* LengthCutoff */ repeat_length_upper_threshold_ && (in_degree.find(from) != in_degree.end()) && (out_degree.find(into) != out_degree.end()) ){

				edge_to_kind_[*e_iter] = LENGTH;
				components.push_back(*e_iter);
			} 
			else{
				singles.push_back(*e_iter);
			}
			
		}
	
		cout << "Number of edges identified by paired info: " << byPairedInfo << endl;
		cout << "Number of edges identified by topology: " << byTopology << endl;
		cout << "Number of edges identified by topology confirmed by paired info: " << confirmedByPairedInfo << endl;
		
		cout << "SINGLES: ";
		for (auto sit = singles.begin(); sit != singles.end(); ++sit){
			cout << gp_.g.int_id(*sit) << ", ";
		}
		cout << endl;

		cout << "COMPONENTS: ";
		for (auto cit = components.begin(); cit != components.end(); ++cit){
			cout << gp_.g.int_id(*cit) << ", ";
		}
		cout << endl; 
		



	}




	template <typename T1, typename T2>
	struct CompareSecond {
		typedef pair<T1, T2> type;
		bool operator ()(type const& a, type const& b) const {
			return a.second < b.second;
	    	}
	};


	bool IfSelfIntersection( const set<EdgeId>& incoming_edges, const set<EdgeId>& outgoing_edges ) {

		vector<EdgeId> v(min(incoming_edges.size(), outgoing_edges.size()));
		vector<EdgeId>::iterator it = set_intersection(incoming_edges.begin(), incoming_edges.end(), outgoing_edges.begin(), outgoing_edges.end(), v.begin());
		                         
		v.resize(it - v.begin());
		if ( v.size() > 0 ) {
			return true;
		}

		return false;

	}

	bool IfCycledComponent( const EdgeId& edge_in, const vector<EdgeId>& component, set<VertexId> visited_vertices ){


		if (visited_vertices.size() == component.size() + 1) return false;
	
		for ( auto edge = component.begin(); edge != component.end(); ++edge ){
				
			if ( gp_.g.EdgeEnd(edge_in) == gp_.g.EdgeStart(*edge) ) {
				if (visited_vertices.find( gp_.g.EdgeEnd(edge_in) ) != visited_vertices.end()) return true;
				return IfCycledComponent( *edge, component, visited_vertices );
				
			}


		}
	
	}		


	void CountDistance( const EdgeId& edge_in, const EdgeId& edge_out, const vector<EdgeId>& component, int& distance ){
	// gets a repetitive component and calculates the length of the longest path in it

		if ( gp_.g.EdgeEnd(edge_in) == gp_.g.EdgeStart(edge_out) ) return;
	
		for ( auto edge = component.begin(); edge != component.end(); ++edge ){
				
				if ( gp_.g.EdgeEnd(edge_in) == gp_.g.EdgeStart(*edge) ) {
				distance += gp_.g.length(*edge);
				CountDistance( *edge, edge_out, component, distance);
				return;
			}
			
		}

	}

	


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

		printf("\n\n");
		sort(probabilities.begin(), probabilities.end());

		printf("%4.10f 0\n", probabilities[0]);
		double pred(probabilities[0]);
		unsigned i(1);
		for ( auto it = probabilities.begin() + 1; it != probabilities.end(); ++it, ++i ) {
			//cout << fabs(*it - pred) << endl;
			if ( fabs(*it - pred) > 0.00000000001 ) {
				//fprintf(file, "%5.2f %5.2f\n", ac, (double) i / (double) probabilities.size());
				printf("%4.10f %4.10f\n", *it, (double) i / (double) probabilities.size());
			}
			pred = *it;
		}
		//fprintf(file,"\n");
		printf("\n\n");
	}



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
						double probability = bm.GetProbablityFromBucketToBucketForDistance (in_bucket, out_bucket, distance, shift) ;
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

/*
	void findClosest(vector<pair<EdgeId, coverage_value>>& incomingEdgesCoverage,
			vector<pair<EdgeId, coverage_value>>& outgoingEdgesCoverage,
			vector<pair<EdgeId,EdgeId>>& pairsOfEdges){

		int Length = min(incomingEdgesCoverage.size(),outgoingEdgesCoverage.size());

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
*/

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

	template <class DetailedCoverage>
	void TraverseComponents( const vector<EdgeId>& components, 
				const vector<EdgeId>& singles, DetailedCoverage& coverage, 
				//set<EdgeId>& usedEdges,
				size_t insert_size,
				vector< vector<EdgeId>> & resolved_paths,
				EdgeQuality<typename graph_pack::graph_t>& quality_labeler) {

		set<EdgeId> visited_edges;
		INFO("Traversing components");
		//FILE* file = fopen("/home/ksenia/path_resolved.log", "w");
		FILE* file = fopen("/home/ksenia/probabilities_22.log", "w");
		
		int allLength = 0;
		for (auto it = components.begin(); it != components.end(); ++it ){
			allLength += gp_.g.length(*it);
		}

		int numberOfComponents = 0;
		int numberOfLargeComponents = 0;
		int numberOfComponentWithDifferentInOutDegree = 0;
		//int numberOfEdgesDetectedByPairedInfoInResolvedComps = 0;

		int filteredByThresholds(0), resolvedPathsNum(0); 
		//int resolvedComponentsByTopology(0), resolvedComponentsByLength(0);//, resolvedComponentsByLengthAndTopology(0);
		vector<int> pathSizes(21,0);
		
        	int number_of_buckets = 20;
		int K_ = cfg::get().K + 1;
		auto bm = BucketMapper<conj_graph_pack::graph_t>(gp_.g, kmer_index_, K_, number_of_buckets);
		bm.InitBuckets( );

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

			vector< vector <double> > transition_probabilities ;
			
			GetComponentInfo(path, incoming_edges, outgoing_edges, coverage, transition_probabilities, quality_labeler );
			
			if (if_loop) {
				if (incoming_edges.size() == 1 && incoming_edges.size() == outgoing_edges.size() )
					pure_tandem += 1;
				else repetitive_tandem += 1;
				cout << "loop!" << endl;
				continue;

			}

			if (path.size() == 0 ) continue;

			if ( incoming_edges.size() == 0 || outgoing_edges.size() == 0) continue;

			ordinal_repeat += 1;

			int longestPathLen = GetLongestPathLength(path);

//			bool containsNotGenomicEdges = false;
			for (auto iter = path.begin(); iter != path.end(); ++iter) {
				if (quality_labeler.quality(*iter) < 0.5) {
				//containsNotGenomicEdges = true;
				}
			}

			if ( ContainsOnlyShortEdges(path) ) {
				INFO("contains only short edges");
				cout << "component: ";
				for ( auto iter = path.begin(); iter != path.end(); ++iter ) {
					cout << gp_.g.int_id(*iter)  << " ";
				}
				cout << endl;


				continue;
			}

			if ( insert_size < (size_t)longestPathLen ) numberOfLargeComponents += 1;

			numberOfComponents += 1;

			if (incoming_edges.size() != outgoing_edges.size() /*&& insert_size < (size_t)longestPathLen */){

				numberOfComponentWithDifferentInOutDegree += 1;
				
				continue;

			}
				

			
	/*		fprintf(file, "resolved path \n");
			for ( auto iter = path.begin(); iter != path.end(); ++iter ) {
				fprintf(file, "%lu ", gp_.g.int_id(*iter));
				fprintf(file, " ");
			}
			fprintf(file, "\nincoming edges: ");
			for ( auto e = incoming_edges.begin(); e != incoming_edges.end(); ++e) {
	
				fprintf(file, "%lu (%5.2f) ", gp_.g.int_id(*e), coverage.GetOutCov(*e));
			}
			fprintf(file,"\n");
			fprintf(file,"outgoing edges: ");
			for ( auto e = outgoing_edges.begin(); e != outgoing_edges.end(); ++e) {
	
				fprintf(file, "%lu (%5.2f) ", gp_.g.int_id(*e), coverage.GetInCov(*e));
			}
			fprintf(file,"\n");

			fprintf(file,"\n");
	*/		
			
			 
/*			vector<pair<EdgeId, coverage_value>> incomingEdgesCoverage, outgoingEdgesCoverage;

			for ( auto inEdge = incoming_edges.begin(); inEdge != incoming_edges.end(); ++inEdge) {
				incomingEdgesCoverage.push_back(make_pair(*inEdge,coverage.GetOutCov(*inEdge)));
			}

			for ( auto outEdge = outgoingEdges.begin(); outEdge != outgoingEdges.end(); ++outEdge) {
				outgoingEdgesCoverage.push_back(make_pair(*outEdge,coverage.GetInCov(*outEdge)));
			}
	
			sort(incomingEdgesCoverage.begin(), incomingEdgesCoverage.end(), CompareSecond<EdgeId, coverage_value>());
			sort(outgoingEdgesCoverage.begin(), outgoingEdgesCoverage.end(), CompareSecond<EdgeId, coverage_value>());
*/			
		/*	INFO("incoming edges");
			for ( auto e = incomingEdgesCoverage.begin(); e != incomingEdgesCoverage.end(); ++e) {
	
				cout << e->first << " " << e->second << endl;
			}
			INFO("outgoing edges");
			for ( auto e = outgoingEdgesCoverage.begin(); e != outgoingEdgesCoverage.end(); ++e) {
	
				cout << e->first << " " << e->second << endl;
			}
		*/
			vector<pair<EdgeId,EdgeId>> pairs_of_edges;


			MatchPairs( incoming_edges, outgoing_edges, pairs_of_edges, path, bm, coverage, quality_labeler, file);

			if ( insert_size < (size_t)longestPathLen )
				if (pairs_of_edges.size() == 0) 
					filteredByThresholds += 1;

			for ( auto edgePair = pairs_of_edges.begin(); edgePair != pairs_of_edges.end(); ++edgePair ){
				path_extend::BidirectionalPath* resolved_path = new path_extend::BidirectionalPath( gp_.g );
				ResolveRepeat( *edgePair, path, *resolved_path );
				/*if (resolved_path->Size() < 3 || resolved_path->Front() == resolved_path->Back() || coverage.GetOutCov(edgePair->first) < threshold_global_ 
					|| coverage.GetInCov(edgePair->second) < threshold_global_) {
					continue;
				}*/

				vector<EdgeId> tempPath = resolved_path->ToVector();
			
				if ( insert_size < (size_t)longestPathLen )
					resolvedPathsNum += 1;
				resolved_paths.push_back( tempPath );
			}

		}
		
		cout << "pure tandems: " << pure_tandem << endl << "repeats + tandems: " << repetitive_tandem << endl << "ordinal repeats: " << ordinal_repeat << endl
		<< "comps with in degree not_equal out degree: " << numberOfComponentWithDifferentInOutDegree << endl;
		//cout << "Number of components : " << numberOfComponents << endl;
		//cout << "Number of components with length exceeding insert size: " << numberOfLargeComponents << endl;
		//cout << "Number of components with length exceeding insert size with different in and out degree: " << numberOfComponentWithDifferentInOutDegree << endl;
		//cout << "Number of components with length exceeding insert size filtered by thresholds: " << filteredByThresholds << endl;
		//cout << "Number of resolved components with length exceeding insert size: " << resolvedPathsNum << endl;
		//cout << "Number of resolved components detected by topology : " << resolvedComponentsByTopology << endl;
		//cout << "Number of resolved components detected by length : " << resolvedComponentsByLength << endl;
		//cout << "Number of edges detected by paired info in resolved comps: " << numberOfEdgesDetectedByPairedInfoInResolvedComps << endl;
		fclose(file);

	}

	template <class EdgesPositionHandlerT> 
	bool match( const path_extend::BidirectionalPath &path, const int currentId, const int currentStart,
	                        EdgesPositionHandlerT& ref_pos) {
		
		//cout << "id: " << path.Size() << " " << currentId << endl;
		if (path.Size() == currentId ){
			return true;
		}
		EdgeId edge = path.At(currentId);
		auto pos_it = ref_pos.edges_positions().find(edge);
		//VERIFY(pos_it != ref_pos.edges_positions().end());

		if ( currentId == path.Size() - 1 && pos_it->second.size() > 1 ){
			//INFO("last is repeat - fail");
			return false;
		}

		bool matched = false;
		for (size_t i = 0; i < pos_it->second.size(); ++i) {
			
			auto start = pos_it->second[i].start();

		      /*cout << "    " << pos_it->second[i].contigId_ << " "
		    			<< pos_it->second[i].start() << " - "
				         << pos_it->second[i].end() << endl;
			*/
			if ( abs(start - currentStart) < 2 ) {

				auto end = pos_it->second[i].end();
				matched = match( path, currentId + 1, end + 1, ref_pos);
			}
		}
		
		return matched;
	}


	bool matchReference( const path_extend::BidirectionalPath& path) {

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

	void dfs( const VertexId& vStart, const VertexId& vEnd, const vector<EdgeId>& component,
		set<EdgeId>& visited,  path_extend::BidirectionalPath& path) {

		for ( auto edge = component.begin(); edge != component.end(); ++edge ){

			if (visited.find(*edge) != visited.end()){
				//cout << "visited contains " << gp_.g.int_id(*edge) << endl;
				continue;
			}
			//cout << gp_.g.int_id(*edge) << endl;
			/*if ( gp_.g.int_id(vStart) == 6402783) {
				cout << "11210071 is start vertex!" << endl;
			}
			if ( gp_.g.int_id(vEnd) == 6402783) {
				cout << "11210071 is end vertex!" << endl;
			}*/
	
			if ( vStart == gp_.g.EdgeStart(*edge) ){
				if ( vEnd == gp_.g.EdgeEnd(*edge) ){
					path.PushBack(*edge);
					visited.insert(*edge);
					return;
				}

				visited.insert(*edge);
				dfs( gp_.g.EdgeEnd(*edge), vEnd, component, visited, path );
				//if (path.Size() > 0){
				path.PushFront(*edge);
				return;
				//}
			}
		}

		return;
	}

	void ResolveRepeat( const pair<EdgeId,EdgeId>& pairOfEdges,
			const vector<EdgeId>& component, path_extend::BidirectionalPath& path ) {
	
		EdgeId incomingEdge = pairOfEdges.first;
		EdgeId outgoingEdge = pairOfEdges.second;

		VertexId vertexStartRepeat = gp_.g.EdgeEnd(incomingEdge); 
		VertexId vertexEndRepeat = gp_.g.EdgeStart(outgoingEdge); 
	
		set<EdgeId> visited;
		dfs(vertexStartRepeat, vertexEndRepeat, component, visited, path);

		path.PushFront(incomingEdge);
		path.PushBack(outgoingEdge);

	
		// here we filter f.e. paths-loops which can be collected when performing dfs because there components go through 
		// vertices which are assigned with edges-loops so that they can be added in the dfs queue 
		if ( path.Size() < 3 ){
			return;
		}
	

	}

	void visit( const EdgeId& edge, set<EdgeId>& visited_edges, set<VertexId>& grey_vertices, 
		vector<EdgeId>& path, const vector<EdgeId>& components, 
		const vector<EdgeId>& singles, vector<EdgeId>& incoming_edges, vector<EdgeId>& outgoing_edges, bool& if_loop ) {

		VertexId edgeStartVertex = gp_.g.EdgeStart(edge);
		VertexId edgeEndVertex = gp_.g.EdgeEnd(edge);
		
		if (visited_edges.find(edge) != visited_edges.end() )
			return;
	
		if ( grey_vertices.find(edgeStartVertex) != grey_vertices.end() ) { //&& grey_vertices.find(edgeEndVertex) != grey_vertices.end() ) {
			if_loop = true;
			return;
		}
		
		grey_vertices.insert(edgeStartVertex);
		//component_vertices.insert(edgeEndVertex);

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
				
				//cout << "visiting edge " << gp_.g.int_id(*component_edge) << endl;
				visit(*component_edge, visited_edges, grey_vertices, path, components, singles, incoming_edges, outgoing_edges, if_loop);
			}


		}

		grey_vertices.insert(edgeStartVertex);

	}


	void GetOtherEdges(path_extend::PathContainer& paths, const set<EdgeId>& used_edges){
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


	void WriteResolved( const set<EdgeId>& used_edges, path_extend::PathContainer& resolved_paths, const string &file_name  ){

		path_extend::ContigWriter cw( gp_.g);
		GetOtherEdges( resolved_paths, used_edges );
		cw.writePaths( resolved_paths, file_name );
	}
};
}

#endif
