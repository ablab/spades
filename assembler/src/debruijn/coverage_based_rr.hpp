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

template <class GraphPack>
class CoverageBasedResolution {
	typedef double coverage_value;
	typedef enum { TOPOLOGY, LENGTH, PAIREDINFO } kind_of_repeat; 
	
	GraphPack *gp;
	const DeBruijnEdgeIndex<EdgeId>& kmer_index_;

	std::map<EdgeId, kind_of_repeat> edge_to_kind_;

	const double threshold_one_list_;
	const double threshold_match_;
	const double threshold_global_;
	const double tandem_lower_threshold_;
	const double tandem_upper_threshold_;
	const double repeat_length_upper_threshold_;
	
	public:
	CoverageBasedResolution( GraphPack *gpack_arg, const DeBruijnEdgeIndex<EdgeId>& kmer_index, double threshold_one_list, double threshold_match, 
				double threshold_global, double tandem_lower_threshold, double tandem_upper_threshold, double repeat_length_upper_threshold) :
											kmer_index_(kmer_index),
											threshold_one_list_(threshold_one_list), 
											threshold_match_(threshold_match),
											threshold_global_(threshold_global), 
											tandem_lower_threshold_(tandem_lower_threshold), 
											tandem_upper_threshold_(tandem_upper_threshold),
											repeat_length_upper_threshold_(repeat_length_upper_threshold){
		gp = gpack_arg;
	}
	
	
	//path without conjugate edges
//	std::vector< PathInfo<typename GraphPack::graph_t> > filteredPaths;


	private:
	void JoinPaths( std::vector<std::vector<EdgeId> >& paths, //std::vector< std::vector<EdgeId> >& paths, 
			std::vector< std::vector<EdgeId> >& resolvedLoops,
			std::vector< std::vector<EdgeId> >& all_paths ) {

		std::map< EdgeId, std::vector<EdgeId> > startEdgeToPath;
		std::map< EdgeId, std::vector<EdgeId> > backEdgeToPath;

		std::vector< std::vector<EdgeId> > bothPaths;
		for (auto path = resolvedLoops.begin(); path != resolvedLoops.end(); ++path) {
			bothPaths.push_back(*path);

		}
		for (auto path = paths.begin(); path != paths.end(); ++path) {
			bothPaths.push_back(*path);
		}

		INFO("before map filling");

//TODO: assert that edge is not more than once first and not more than once last in these paths.	
		for (auto path = bothPaths.begin(); path != bothPaths.end(); ++path) {
		//	std::cout << gp->g.int_id(path->front()) << std::endl;
		//	VERIFY(startEdgeToPath.find(path->front()) == startEdgeToPath.end());
			startEdgeToPath[path->front()] = *path;
		}

		for (auto path = bothPaths.begin(); path != bothPaths.end(); ++path) {
		//	VERIFY(backEdgeToPath.find(path->back()) == backEdgeToPath.end());
			backEdgeToPath[path->back()] = *path;
		}

		for (auto path = bothPaths.begin(); path != bothPaths.end(); ++path) {
			if ( backEdgeToPath.find( path->front() ) != backEdgeToPath.end()) {
				continue;
			}

			/*INFO("path before");
			for ( auto e = path->begin(); e != path->end(); ++e ){
				std::cout << gp->g.int_id(*e) <<  " ";
			}
			std::cout << std::endl; */
			bool updated = true;
			while (updated) {
				
				auto foundPath = startEdgeToPath.find( path->back() );
				if (foundPath != startEdgeToPath.end() ) {
					//INFO("found path before");
					/*for ( auto e = foundPath->second.begin(); e != foundPath->second.end(); ++e ){
						std::cout << gp->g.int_id(*e) <<  " ";
					}
					std::cout << std::endl;*/
					path->insert(path->end(), foundPath->second.begin() + 1, foundPath->second.end());
				}

				else {
					updated = false;
				}

			}
			/*INFO("path after");
			for ( auto e = path->begin(); e != path->end(); ++e ){
				std::cout << gp->g.int_id(*e) <<  " ";
			}
			std::cout << std::endl;*/
			all_paths.push_back(*path);
		}
		INFO("out of path joining");
	
	}



	bool VerifyComponent( std::vector<EdgeId>& incoming_edges,
			 std::vector<EdgeId>& outgoing_edges,
		 	std::vector<EdgeId>& component ){

		if ( incoming_edges.size() == outgoing_edges.size() ) {
			return true;
		}

		int diff = incoming_edges.size() - outgoing_edges.size();

		if ( diff < 0 ) {
			
			int counter = 0;
			for ( auto edge = incoming_edges.begin(); edge != incoming_edges.end(); ++edge ) {

				if ( gp->g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
					counter += 1;
				}
			}

			if ( counter == -diff ) {
				
				std::cout << "INCOMING COMPONENT UPDATED" << std::endl;
				/*for ( auto edge = incomingEdges.begin(); edge != incomingEdges.end(); ++edge ) {

					if ( gp->g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
						components.push_back(*edge);	
					}
				}*/
			}
		
		}

		if ( diff > 0 ) {
			
			int counter = 0;
			for ( auto edge = outgoing_edges.begin(); edge != outgoing_edges.end(); ++edge ) {

				if ( gp->g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
					counter += 1;
				}
			}

			if ( counter == -diff ) {
				
				std::cout << "OUTGOING COMPONENT UPDATED" << std::endl;
				/*for ( auto edge = incomingEdges.begin(); edge != incomingEdges.end(); ++edge ) {

					if ( gp->g.length(*edge) <= cfg::get().rr.max_repeat_length ) {
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
					EdgeLabelHandler<typename GraphPack::graph_t>& labels_after,
					EdgeQuality<typename GraphPack::graph_t>& quality_labeler,
					PairedInfoIndexT<typename GraphPack::graph_t> & clustered_index,
					std::vector< PathInfo<typename GraphPack::graph_t> >& filtered_paths
					//std::set<EdgeId>& prohibitedEdges,
					//std::vector< std::vector<EdgeId>>& resolvedLoops ) 
					) {

		auto filter = LoopFilter<GraphPack, DetailedCoverage>(*gp, coverage, tandem_lower_threshold_, tandem_upper_threshold_, repeat_length_upper_threshold_);
		filter.get_loopy_components( quality_labeler );


		INFO("Resolving repeats by coverage...");

		std::cout << "prohibited edges" << std::endl;
		for (auto e_iter = filter.prohibitedEdges.begin(); e_iter!= filter.prohibitedEdges.end(); ++e_iter){
			std::cout << gp->g.int_id(*e_iter) << ", ";
		}
		std::cout << std::endl;

		std::vector<EdgeId> components, singles;
		std::vector<EdgeId> componentsRef, singlesRef;
		INFO("Getting components...");
		GetComponents(  components, singles, labels_after, quality_labeler, clustered_index, filter.prohibitedEdges );
		GetComponentsWithReference(  componentsRef, singlesRef, quality_labeler, filter.prohibitedEdges );

		
		int fp(0), fn(0);

		std::cout << "in components (size: " << components.size() << ") but not in componentsRef: " << std::endl;
		for (auto it = components.begin(); it != components.end(); ++it) {

			if (std::find(componentsRef.begin(), componentsRef.end(), *it) == componentsRef.end() 
				&& (quality_labeler.quality(*it) > 0.5) ) {
				fp += 1;
				std::cout << gp->g.int_id(*it);
				if ( edge_to_kind_[*it] == TOPOLOGY ) {
					  std::cout << " (TOPOLOGY) , ";
				}
				if ( edge_to_kind_[*it] == LENGTH ) {
					std::cout << " (LENGTH) , ";
				}

			}

		}
		std::cout << std::endl;
		
	
		std::cout << "in componentsRef (size: " << componentsRef.size() << ") but not in components: " << std::endl;
		for (auto it = componentsRef.begin(); it != componentsRef.end(); ++it) {

			if (std::find(components.begin(), components.end(), *it) == components.end()) {
				fn += 1;
				std::cout << gp->g.int_id(*it) << ", ";
			}

		}
		std::cout << std::endl;
		std::cout << "False positives: " << (double) fp / components.size() << "False negatives: " << (double) fn / (fn + singles.size()) << std::endl;
			

		//path with conjugate edges
		std::vector< std::vector<EdgeId> > resolved_paths;
	
		//getComponents( gp, components, singles, quality_labeler, unresolvedLoops );

		INFO("Traversing graph...");
		TraverseComponents( components, singles, coverage, insert_size, resolved_paths, quality_labeler );

		std::set<EdgeId> used_edges;

		std::cout << "Paths before joining: " << std::endl;
		for ( auto p = resolved_paths.begin(); p != resolved_paths.end(); ++p) {
			//fprintf(file, "resolved path \n");
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				std::cout << gp->g.int_id(*iter) << " (" << gp->g.int_id(gp->g.EdgeStart(*iter) ) << "," << gp->g.int_id(gp->g.EdgeEnd(*iter) ) << ") ";
				//fprintf(file, "%d ", gp->g.int_id(*iter));
				//fprintf(file, " ");
			}
			std::cout << std::endl;
			//fprintf(file,"\n");
		}
		std::cout << "Loops before joining: " << std::endl;
		for ( auto p = filter.resolvedLoops.begin(); p != filter.resolvedLoops.end(); ++p) {
			//fprintf(file, "resolved path \n");
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				std::cout << gp->g.int_id(*iter) << "( " << gp->g.length(*iter) << ") ";
				//fprintf(file, "%d ", gp->g.int_id(*iter));
				//fprintf(file, " ");
			}
			std::cout << std::endl;
			//fprintf(file,"\n");
		}

		std::vector< std::vector<EdgeId> > all_paths;
		JoinPaths(resolved_paths, filter.resolvedLoops, all_paths);
		FilterConjugate( used_edges, all_paths, filtered_paths);
		//std::cout << "before filtering size " << allPaths.size() << " filtered size: " << filteredPaths.size() << std::endl;

		
		for ( auto p = all_paths.begin(); p != all_paths.end(); ++p) {
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				std::cout << gp->g.int_id(*iter) << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "-------------------------" << std::endl;
		for ( auto p = filtered_paths.begin(); p != filtered_paths.end(); ++p) {
			//fprintf(file, "resolved path \n");
			for ( auto iter = p->getPath().begin(); iter != p->getPath().end(); ++iter ) {
				std::cout << gp->g.int_id(*iter) << "( " << gp->g.length(*iter) << "; " << coverage.GetInCov(*iter) << " " << coverage.GetOutCov(*iter) << ") ";
				//fprintf(file, "%d ", gp->g.int_id(*iter));
				//fprintf(file, " ");
			}
			std::cout << std::endl;
			//fprintf(file,"\n");
		}

		//fclose(file);
		std::string file_name = cfg::get().output_dir + "resolved_by_coverage.fasta";
		//INFO("Writing result");

		path_extend::PathContainer paths_output;
		for ( auto p = filtered_paths.begin(); p != filtered_paths.end(); ++p) {
			path_extend::BidirectionalPath* bidirectional_path = new path_extend::BidirectionalPath( gp->g );
			path_extend::BidirectionalPath* conjugate_path = new path_extend::BidirectionalPath( gp->g );
			auto tmpPath = p->getPath();
			for (auto it = tmpPath.begin(); it != tmpPath.end(); ++it ){
				
					bidirectional_path->PushBack(*it);
					EdgeId cedge = gp->g.conjugate(*it);
					conjugate_path->PushFront( cedge );
			} 

			paths_output.AddPair( bidirectional_path, conjugate_path );
		}

		WriteResolved( used_edges, paths_output, file_name);
	}

	private:

	void FilterConjugate( std::set<EdgeId>& used_edges,
				const std::vector< std::vector<EdgeId> > & paths,
				std::vector< PathInfo<typename GraphPack::graph_t> >& filtered_paths) {


		INFO("filtering conjugate edges");
		for ( auto path = paths.begin(); path != paths.end(); ++path) {

			bool ifInsert = true;
			for (auto e = path->begin(); e != path->end(); ++e) {
				if ( gp->g.conjugate(*e) == *e ) continue; 
				if ( used_edges.find(gp->g.conjugate(*e)) != used_edges.end() ) {
//TODO:: this is not true, if we have autoreverse edge.					
					ifInsert = false;
					break;
				}

			}
			if (ifInsert) {
				/*if (! gp->edge_pos.IsConsistentWithGenome(*path)) {
					for (auto iter = path->begin(); iter != path->end(); ++iter) {
						auto positions = gp->edge_pos.GetEdgePositions(*iter);
				}*/

				filtered_paths.push_back(PathInfo<typename GraphPack::graph_t>(*path));
				for (auto e = path->begin(); e != path->end(); ++e) 
					used_edges.insert(*e);
				}
		}
		
		//INFO("inserting paths into set of the used edges");
		for ( auto path = filtered_paths.begin(); path != filtered_paths.end(); ++path ) {
			auto p = path->getPath();
			//std::cout << "size of p: " << p.size() << std::endl;
			for (auto e = p.begin(); e != p.end(); ++e) {
				used_edges.insert(*e);
			}
		}

		//INFO("out of filtering");
	}

	void GetComponentsWithReference( std::vector<EdgeId>& components, std::vector<EdgeId>& singles,
					EdgeQuality<typename GraphPack::graph_t>& quality_labeler,
					std::set<EdgeId>& prohibitedEdges ){

		INFO("Finding Components With Paired Info");
		for (auto iter = gp->g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {

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

		auto edgeStart = gp->g.EdgeStart(edge);
		auto edgeEnd = gp->g.EdgeEnd(edge);
		auto outgoingFromStart = gp->g.OutgoingEdges(edgeStart);
		//auto incomingToEnd = gp->g.IncomingEdges(edgeEnd);
 
		for (auto e = outgoingFromStart.begin(); e != outgoingFromStart.end(); ++e){
			if (*e == edge) continue;
			if (gp->g.EdgeEnd(*e) == edgeEnd) return true;

		}

		return false;
	}

	template< class Graph>
	bool CheckIfComponentByPairedInfo( EdgeId edge, PairedInfoIndexT<Graph>& clustered_index, std::set<EdgeId>& prohibitedEdges ) {

		io::SequencingLibrary<debruijn_config::DataSetData> lib;
		auto improver = PairInfoImprover<Graph>(gp->g, clustered_index,lib);
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
					//std::cout << "Inconsistent for " << gp->g.int_id(edge) << ": " << gp->g.int_id(e1) << " " << gp->g.int_id(e2) << std::endl;
					return true;
				}
			}
		}

		return false;
	}

	void GetComponents( std::vector<EdgeId>& components, std::vector<EdgeId>& singles, 
				EdgeLabelHandler<typename GraphPack::graph_t>& labels_after,
				EdgeQuality<Graph>& quality_labeler,
				PairedInfoIndexT<Graph>& clustered_index,
				std::set<EdgeId>& prohibitedEdges ){

		typedef int times;
		std::map<VertexId, times> out_degree;
		std::map<VertexId, times> in_degree;
	/*	double LengthCutoff = 0.0;

		int numEdges = 0;
		for (auto e_iter = gp->g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
			LengthCutoff += gp->g.length(*e_iter);	
			numEdges += 1;
		}

		LengthCutoff = 0.25 * LengthCutoff / numEdges;	
	*/

//		std::cout << "Length Cutoff: "  << LengthCutoff << std::endl;
		INFO("Getting Components");
		INFO("Counting degrees of vertices");
		for (auto e_iter = gp->g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
		
			//if ( prohibitedEdges.find(*e_iter) != prohibitedEdges.end() ) continue;
			VertexId from = gp->g.EdgeStart(*e_iter);
			VertexId into = gp->g.EdgeEnd(*e_iter);

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
		for (auto e_iter = gp->g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
			LengthCutoff += gp->g.length(*e_iter);	
			numEdges += 1;
		}

		std::cout << "Length cutoff: " << 0.25 * LengthCutoff / numEdges << std::endl;
		*/
		int byPairedInfo = 0;
		int byTopology = 0;
		int confirmedByPairedInfo = 0;
		for (auto e_iter = gp->g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){

			if ( prohibitedEdges.find(*e_iter) != prohibitedEdges.end() ){ 
				continue;
			}

			VertexId from = gp->g.EdgeStart(*e_iter);
			VertexId into = gp->g.EdgeEnd(*e_iter);

			if ( gp->g.length(*e_iter) >= cfg::get().rr.max_repeat_length ) {

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
				//	std::cout << "Component Edge Detected By Paired Info: " << gp->g.int_id(*e_iter) << std::endl;
				byPairedInfo += 1;
				edge_to_kind_[*e_iter] = PAIREDINFO;
			}*/

			else if( gp->g.length(*e_iter) < /* LengthCutoff */ repeat_length_upper_threshold_ && (in_degree.find(from) != in_degree.end()) && (out_degree.find(into) != out_degree.end()) ){

				edge_to_kind_[*e_iter] = LENGTH;
				components.push_back(*e_iter);
			} 
			else{
				singles.push_back(*e_iter);
			}
			
		}
	
		std::cout << "Number of edges identified by paired info: " << byPairedInfo << std::endl;
		std::cout << "Number of edges identified by topology: " << byTopology << std::endl;
		std::cout << "Number of edges identified by topology confirmed by paired info: " << confirmedByPairedInfo << std::endl;
		
		std::cout << "SINGLES: ";
		for (auto sit = singles.begin(); sit != singles.end(); ++sit){
			std::cout << gp->g.int_id(*sit) << ", ";
		}
		std::cout << std::endl;

		std::cout << "COMPONENTS: ";
		for (auto cit = components.begin(); cit != components.end(); ++cit){
			std::cout << gp->g.int_id(*cit) << ", ";
		}
		std::cout << std::endl; 
		



	}




	template <typename T1, typename T2>
	struct CompareSecond {
		typedef pair<T1, T2> type;
		bool operator ()(type const& a, type const& b) const {
			return a.second < b.second;
	    	}
	};


	bool IfSelfIntersection( const std::set<EdgeId>& incoming_edges, const std::set<EdgeId>& outgoing_edges ) {

		std::vector<EdgeId> v(min(incoming_edges.size(), outgoing_edges.size()));
		std::vector<EdgeId>::iterator it = std::set_intersection(incoming_edges.begin(), incoming_edges.end(), outgoing_edges.begin(), outgoing_edges.end(), v.begin());
		                         
		v.resize(it - v.begin());
		if ( v.size() > 0 ) {
			return true;
		}

		return false;

	}

	bool IfCycledComponent( const EdgeId& edge_in, const std::vector<EdgeId>& component, std::set<VertexId> visited_vertices ){


		if (visited_vertices.size() == component.size() + 1) return false;
	
		for ( auto edge = component.begin(); edge != component.end(); ++edge ){
				
			if ( gp->g.EdgeEnd(edge_in) == gp->g.EdgeStart(*edge) ) {
				if (visited_vertices.find( gp->g.EdgeEnd(edge_in) ) != visited_vertices.end()) return true;
				return IfCycledComponent( *edge, component, visited_vertices );
				
			}


		}
	
	}		


	void CountDistance( const EdgeId& edge_in, const EdgeId& edge_out, const std::vector<EdgeId>& component, int& distance ){
	// gets a repetitive component and calculates the length of the longest path in it

		if ( gp->g.EdgeEnd(edge_in) == gp->g.EdgeStart(edge_out) ) return;
	
		for ( auto edge = component.begin(); edge != component.end(); ++edge ){
				
				if ( gp->g.EdgeEnd(edge_in) == gp->g.EdgeStart(*edge) ) {
				distance += gp->g.length(*edge);
				CountDistance( *edge, edge_out, component, distance);
				return;
			}
			
		}

	}

	template <class DetailedCoverage>
	bool MatchPairs( std::vector<EdgeId>& incoming_edges,
			 std::vector<EdgeId>& outgoing_edges,
			 std::vector<std::pair<EdgeId,EdgeId>>& pairs_of_edges,
			 const std::vector<EdgeId>& component,
			 BucketMapper<Graph> &bm,
			 DetailedCoverage& coverage,
			 EdgeQuality<typename GraphPack::graph_t>& quality_labeler)  {


//		if (incomingEdges.size() > 5 || outgoingEdges.size() > 5) return false;
		double shift = 25;

		std::vector< std::vector <double> > transition_probabilities ;
		for ( unsigned i = 0; i < incoming_edges.size(); i++) {
			 transition_probabilities.push_back(std::vector<double>(outgoing_edges.size(),1));
		}

		int in_edge_counter(0);


		//std::cout << "component size: " << component.size() << std::endl;
		for ( auto in_edge = incoming_edges.begin(); in_edge != incoming_edges.end(); ++in_edge, ++in_edge_counter ) {
	
			//std::cout << "incoming edge " << gp->g.int_id(*in_edge) << std::endl;
	
			double in_cov = coverage.GetOutCov(*in_edge);
			int in_bucket = bm.GetCoverageBucket(in_cov);
			 
			 	int out_edge_counter(0);
				for ( auto out_edge = outgoing_edges.begin(); out_edge != outgoing_edges.end(); ++out_edge, ++out_edge_counter ) {
				
						
						//std::cout << "outgoing edge " << gp->g.int_id(*out_edge) << std::endl;
						double out_cov = coverage.GetInCov(*out_edge);
						//std::cout << "out_cov: " << out_cov << std::endl;
						int out_bucket = bm.GetCoverageBucket(out_cov);
						//std::cout << "out_bucket: " << out_bucket << std::endl;

						int distance(0);
						CountDistance(*in_edge, *out_edge, component, distance);

						//std::cout << distance << " " << in_bucket << " " << out_bucket << std::endl;
						double probability = bm.GetProbablityFromBucketToBucketForDistance (in_bucket, out_bucket, distance, shift) ;
						//std::cout << probability << std::endl;
						transition_probabilities[in_edge_counter][out_edge_counter] = probability;
		        	} 
				//std::cout << std::endl;
		}

		GetComponentInfo(component, incoming_edges, outgoing_edges, coverage, transition_probabilities, quality_labeler );

			
		unsigned k = 0;
		std::set<unsigned> matched;
		std::cout << "MATCH: " << std::endl;
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
				if ( std::abs(max_prob - (*vec)[i]) < 0.1 )
					if_matches = false;
			}

			if (if_matches)
			for (unsigned j = 0; j < vec->size(); ++j) {
				if (j == k) continue;
				if ( std::abs(max_prob - transition_probabilities[j][max_index]) < 0.1 )
					if_matches = false;
			}

			if (if_matches) {
				std::cout << k << " " <<  max_index << std::endl;
				matched.insert(k);
				matched.insert(max_index);
				pairs_of_edges.push_back(std::make_pair(incoming_edges[k], outgoing_edges[max_index]) );
			}

		}
		
		if (pairs_of_edges.size() > 0) return true;

		return false;
		
	}
/*
	void findClosest(std::vector<std::pair<EdgeId, coverage_value>>& incomingEdgesCoverage,
			std::vector<std::pair<EdgeId, coverage_value>>& outgoingEdgesCoverage,
			std::vector<std::pair<EdgeId,EdgeId>>& pairsOfEdges){

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
			pairsOfEdges.push_back( std::make_pair(incomingEdgesCoverage[i].first, outgoingEdgesCoverage[i].first)  );
		}

	}
*/

	bool ContainsSmallLoop( const std::vector<EdgeId>& path){
		
		if ( path.size() == 1 ) return true;

		if ( path.size() == 2 ) {

			if (gp->g.EdgeStart(path[0]) == gp->g.EdgeEnd(path[1]) && gp->g.EdgeStart(path[1]) == gp->g.EdgeEnd(path[0]))
				return true;
		}

		return false;
	}

	bool ContainsOnlyShortEdges( const std::vector<EdgeId>& path){

		for ( auto it = path.begin(); it != path.end(); ++it ) {
			
			if (gp->g.length(*it) >= repeat_length_upper_threshold_ || edge_to_kind_[*it] == TOPOLOGY )
				return false;
			
		}
		return true;
	}


	void bfs ( const EdgeId& edge,  std::set<EdgeId>& visited_edges, const std::vector<EdgeId>& component, int& curLen, int& maxPathLen) {

		visited_edges.insert(edge);
		auto incoming_edges = gp->g.IncomingEdges(gp->g.EdgeStart(edge));

		for ( auto e = incoming_edges.begin(); e != incoming_edges.end(); ++e) {
									
		
			if ( std::find(component.begin(), component.end(), *e) != component.end() && visited_edges.find(*e) == visited_edges.end() ){
				curLen += gp->g.length(*e);
				if (curLen > maxPathLen) maxPathLen = curLen;
				bfs(*e, visited_edges, component, curLen, maxPathLen);
			}

		}

		auto outgoingEdges = gp->g.OutgoingEdges(gp->g.EdgeEnd(edge));
		for ( auto e = outgoingEdges.begin(); e != outgoingEdges.end(); ++e) {
	
			if ( std::find(component.begin(), component.end(), *e) != component.end() && visited_edges.find(*e) == visited_edges.end() ){
					
					curLen += gp->g.length(*e);
					if (curLen > maxPathLen) maxPathLen = curLen;
					bfs(*e, visited_edges, component, curLen, maxPathLen);
			}
		}

	}

	int GetLongestPathLength( const std::vector<EdgeId>& component ){
	// gets a repetitive component and calculates the length of the longest path in it

		std::set<EdgeId> visited_edges;
		std::vector<std::vector<EdgeId>> paths;

		int maxPathLen = gp->g.length(component[0]);
		for ( auto edge = component.begin(); edge != component.end(); ++edge ){

			if (visited_edges.find(*edge) != visited_edges.end()) continue;
			int curLen = gp->g.length(*edge);
			visited_edges.insert(*edge);
			bfs(*edge, visited_edges, component, curLen, maxPathLen);	
			
		}

		return maxPathLen;

	}


	bool CheckRepeatDetection( const std::vector<EdgeId>& component, 
				const std::vector<EdgeId>& incoming_edges, 
				const std::vector<EdgeId>& outgoing_edges,
				const EdgeQuality<typename GraphPack::graph_t>& quality_labeler ) {

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
	void GetComponentInfo(const std::vector<EdgeId>& component, const std::vector<EdgeId>& incoming_edges, const std::vector<EdgeId>& outgoing_edges, const DetailedCoverage& coverage, 
		const std::vector< std::vector <double> >& transition_probabilities, const EdgeQuality<typename GraphPack::graph_t>& quality_labeler ) {

			std::cout << "Component: " << std::endl;
			for ( auto iter = component.begin(); iter != component.end(); ++iter ) {
				std::cout << gp->g.int_id(*iter)  << " edge length: " << gp->g.length(*iter) << 
							" average edge coverage " << gp->g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << " "; 
				
				if (edge_to_kind_[*iter] == TOPOLOGY) 
					std::cout << "TOPOLOGY" << std::endl;
				else if  (edge_to_kind_[*iter] == LENGTH )
					std::cout << "LENGTH" << std::endl;
				else if (edge_to_kind_[*iter] == PAIREDINFO )
					std::cout << "PAIREDINFO" << std::endl;
		}
			std::cout << std::endl;
			std::cout << "incoming edges: " << std::endl;
			for ( auto iter = incoming_edges.begin(); iter != incoming_edges.end(); ++iter ) {
				std::cout << gp->g.int_id(*iter)  << " edge length: " << gp->g.length(*iter) << " outgoing edge coverage: " << coverage.GetOutCov(*iter) << 
							" average edge coverage " << gp->g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << std::endl;
			}
			std::cout << std::endl;
			std::cout << "outgoing edges: " << std::endl;
			for ( auto iter = outgoing_edges.begin(); iter != outgoing_edges.end(); ++iter ) {
				std::cout << gp->g.int_id(*iter)  << " edge length: " << gp->g.length(*iter) << " incoming edge coverage: " << coverage.GetInCov(*iter) << 
							" average edge coverage " << gp->g.coverage(*iter) << " quality: " << quality_labeler.quality(*iter) << std::endl;
			}
			std::cout << std::endl;
			
			bool correct_component = checkRepeatDetection( component, incoming_edges, outgoing_edges, quality_labeler  );

			if (!correct_component) {
				std::cout << "repeat is detected incorrectly" << std::endl;
			}
			if (transition_probabilities.size() > 0) {
				std::cout << "transition probabilities" << incoming_edges.size() << "x" << outgoing_edges.size() << ":" << std::endl;
				for (auto vec = transition_probabilities.begin(); vec != transition_probabilities.end(); ++vec ) {
			
					for ( auto prob = vec->begin(); prob != vec->end(); ++prob ) {

						printf("%4.2f ", *prob);
					}
					std::cout << std::endl;
				}
			}
			std::cout << std::endl;

	}

	template <class DetailedCoverage>
	void TraverseComponents( const std::vector<EdgeId>& components, 
				const std::vector<EdgeId>& singles, DetailedCoverage& coverage, 
				//std::set<EdgeId>& usedEdges,
				size_t insert_size,
				std::vector< std::vector<EdgeId>> & resolved_paths,
				EdgeQuality<typename GraphPack::graph_t>& quality_labeler) {

		std::set<EdgeId> visited_edges;
		INFO("Traversing components");
		//FILE* file = fopen("/home/ksenia/path_resolved.log", "w");
		
		int allLength = 0;
		for (auto it = components.begin(); it != components.end(); ++it ){
			allLength += gp->g.length(*it);
		}

		int numberOfComponents = 0;
		int numberOfLargeComponents = 0;
		int numberOfComponentWithDifferentInOutDegree = 0;
		//int numberOfEdgesDetectedByPairedInfoInResolvedComps = 0;

		int filteredByThresholds(0), resolvedPathsNum(0); 
		//int resolvedComponentsByTopology(0), resolvedComponentsByLength(0);//, resolvedComponentsByLengthAndTopology(0);
		std::vector<int> pathSizes(21,0);
		
        	int number_of_buckets = 20;
		int K_ = cfg::get().K + 1;
		auto bm = BucketMapper<conj_graph_pack::graph_t>(gp->g, kmer_index_, K_, number_of_buckets);
		bm.InitBuckets( );

		int pure_tandem(0), repetitive_tandem(0), ordinal_repeat(0);
		for ( auto edge = components.begin(); edge != components.end(); ++edge ) {
			
			if ( visited_edges.find(*edge) != visited_edges.end() ){
				continue;
			}


			std::vector<EdgeId> incoming_edges, outgoing_edges;
			std::vector<EdgeId> path;
		
			std::set<VertexId> component_vertices;
			bool if_loop = false;
			visit(*edge, visited_edges, component_vertices, path, components, singles, incoming_edges, outgoing_edges, if_loop);

			std::vector< std::vector <double> > transition_probabilities ;
			
			GetComponentInfo(path, incoming_edges, outgoing_edges, coverage, transition_probabilities, quality_labeler );
			
			if (if_loop) {
				if (incoming_edges.size() == 1 && incoming_edges.size() == outgoing_edges.size() )
					pure_tandem += 1;
				else repetitive_tandem += 1;
				std::cout << "loop!" << std::endl;
				continue;

			}

			if (path.size() == 0 ) continue;

			if ( incoming_edges.size() == 0 || outgoing_edges.size() == 0) continue;

			ordinal_repeat += 1;

			int longestPathLen = getLongestPathLength(path);

			bool containsNotGenomicEdges = false;
			for (auto iter = path.begin(); iter != path.end(); ++iter) {
				if (quality_labeler.quality(*iter) < 0.5) {
					containsNotGenomicEdges = true;
				}
			}

			if ( containsOnlyShortEdges(path) ) {
				INFO("contains only short edges");
				std::cout << "component: ";
				for ( auto iter = path.begin(); iter != path.end(); ++iter ) {
					std::cout << gp->g.int_id(*iter)  << " ";
				}
				std::cout << std::endl;


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
				fprintf(file, "%lu ", gp->g.int_id(*iter));
				fprintf(file, " ");
			}
			fprintf(file, "\nincoming edges: ");
			for ( auto e = incoming_edges.begin(); e != incoming_edges.end(); ++e) {
	
				fprintf(file, "%lu (%5.2f) ", gp->g.int_id(*e), coverage.GetOutCov(*e));
			}
			fprintf(file,"\n");
			fprintf(file,"outgoing edges: ");
			for ( auto e = outgoing_edges.begin(); e != outgoing_edges.end(); ++e) {
	
				fprintf(file, "%lu (%5.2f) ", gp->g.int_id(*e), coverage.GetInCov(*e));
			}
			fprintf(file,"\n");

			fprintf(file,"\n");
	*/		
			
			 
/*			std::vector<std::pair<EdgeId, coverage_value>> incomingEdgesCoverage, outgoingEdgesCoverage;

			for ( auto inEdge = incoming_edges.begin(); inEdge != incoming_edges.end(); ++inEdge) {
				incomingEdgesCoverage.push_back(std::make_pair(*inEdge,coverage.GetOutCov(*inEdge)));
			}

			for ( auto outEdge = outgoingEdges.begin(); outEdge != outgoingEdges.end(); ++outEdge) {
				outgoingEdgesCoverage.push_back(std::make_pair(*outEdge,coverage.GetInCov(*outEdge)));
			}
	
			sort(incomingEdgesCoverage.begin(), incomingEdgesCoverage.end(), CompareSecond<EdgeId, coverage_value>());
			sort(outgoingEdgesCoverage.begin(), outgoingEdgesCoverage.end(), CompareSecond<EdgeId, coverage_value>());
*/			
		/*	INFO("incoming edges");
			for ( auto e = incomingEdgesCoverage.begin(); e != incomingEdgesCoverage.end(); ++e) {
	
				std::cout << e->first << " " << e->second << std::endl;
			}
			INFO("outgoing edges");
			for ( auto e = outgoingEdgesCoverage.begin(); e != outgoingEdgesCoverage.end(); ++e) {
	
				std::cout << e->first << " " << e->second << std::endl;
			}
		*/
			std::vector<std::pair<EdgeId,EdgeId>> pairs_of_edges;


			MatchPairs( incoming_edges, outgoing_edges, pairs_of_edges, path, bm, coverage, quality_labeler);

			if ( insert_size < (size_t)longestPathLen )
				if (pairs_of_edges.size() == 0) 
					filteredByThresholds += 1;

			for ( auto edgePair = pairs_of_edges.begin(); edgePair != pairs_of_edges.end(); ++edgePair ){
				path_extend::BidirectionalPath* resolved_path = new path_extend::BidirectionalPath( gp->g );
				ResolveRepeat( *edgePair, path, *resolved_path );
				/*if (resolved_path->Size() < 3 || resolved_path->Front() == resolved_path->Back() || coverage.GetOutCov(edgePair->first) < threshold_global_ 
					|| coverage.GetInCov(edgePair->second) < threshold_global_) {
					continue;
				}*/

				std::vector<EdgeId> tempPath = resolved_path->ToVector();
			
				if ( insert_size < (size_t)longestPathLen )
					resolvedPathsNum += 1;
				resolved_paths.push_back( tempPath );
			}

		}
		
		std::cout << "pure tandems: " << pure_tandem << std::endl << "repeats + tandems: " << repetitive_tandem << std::endl << "ordinal repeats: " << ordinal_repeat << std::endl
		<< "comps with in degree not_equal out degree: " << numberOfComponentWithDifferentInOutDegree << std::endl;
		//std::cout << "Number of components : " << numberOfComponents << std::endl;
		//std::cout << "Number of components with length exceeding insert size: " << numberOfLargeComponents << std::endl;
		//std::cout << "Number of components with length exceeding insert size with different in and out degree: " << numberOfComponentWithDifferentInOutDegree << std::endl;
		//std::cout << "Number of components with length exceeding insert size filtered by thresholds: " << filteredByThresholds << std::endl;
		//std::cout << "Number of resolved components with length exceeding insert size: " << resolvedPathsNum << std::endl;
		//std::cout << "Number of resolved components detected by topology : " << resolvedComponentsByTopology << std::endl;
		//std::cout << "Number of resolved components detected by length : " << resolvedComponentsByLength << std::endl;
		//std::cout << "Number of edges detected by paired info in resolved comps: " << numberOfEdgesDetectedByPairedInfoInResolvedComps << std::endl;
		fclose(file);

	}

	template <class EdgesPositionHandlerT> 
	bool match( const path_extend::BidirectionalPath &path, const int currentId, const int currentStart,
	                        EdgesPositionHandlerT& ref_pos) {
		
		//std::cout << "id: " << path.Size() << " " << currentId << std::endl;
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

		      /*std::cout << "    " << pos_it->second[i].contigId_ << " "
		    			<< pos_it->second[i].start() << " - "
				         << pos_it->second[i].end() << std::endl;
			*/
			if ( abs(start - currentStart) < 2 ) {

				auto end = pos_it->second[i].end();
				matched = match( path, currentId + 1, end + 1, ref_pos);
			}
		}
		
		return matched;
	}


	bool matchReference( const path_extend::BidirectionalPath& path) {

		auto ref_pos = gp->edge_pos;
		EdgeId edge = path.At(0);
		auto pos_it = ref_pos.edges_positions().find(edge);
		//VERIFY(pos_it != ref_pos.edges_positions().end();

		if ( pos_it->second.size() == 1 ){
		
			auto nextStart = pos_it->second[0].end();
			return match( path, 1, nextStart + 1, ref_pos );
		}
		//INFO("first is repeat - fail or uncovered");
		std::cout << pos_it->second.size() << std::endl;
		for (size_t i = 0; i < pos_it->second.size(); i++) {
		      std::cout << "    " << pos_it->second[i].contigId_ << " "
		    			<< pos_it->second[i].start() << " - "
				         << pos_it->second[i].end() << std::endl;
		}
		return false;
	}

	void dfs( const VertexId& vStart, const VertexId& vEnd, const std::vector<EdgeId>& component,
		std::set<EdgeId>& visited,  path_extend::BidirectionalPath& path) {

		for ( auto edge = component.begin(); edge != component.end(); ++edge ){

			if (visited.find(*edge) != visited.end()){
				//std::cout << "visited contains " << gp->g.int_id(*edge) << std::endl;
				continue;
			}
			//std::cout << gp->g.int_id(*edge) << std::endl;
			/*if ( gp->g.int_id(vStart) == 6402783) {
				std::cout << "11210071 is start vertex!" << std::endl;
			}
			if ( gp->g.int_id(vEnd) == 6402783) {
				std::cout << "11210071 is end vertex!" << std::endl;
			}*/
	
			if ( vStart == gp->g.EdgeStart(*edge) ){
				if ( vEnd == gp->g.EdgeEnd(*edge) ){
					path.PushBack(*edge);
					visited.insert(*edge);
					return;
				}

				visited.insert(*edge);
				dfs( gp->g.EdgeEnd(*edge), vEnd, component, visited, path );
				//if (path.Size() > 0){
				path.PushFront(*edge);
				return;
				//}
			}
		}

		return;
	}

	void ResolveRepeat( const std::pair<EdgeId,EdgeId>& pairOfEdges,
			const std::vector<EdgeId>& component, path_extend::BidirectionalPath& path ) {
	
		EdgeId incomingEdge = pairOfEdges.first;
		EdgeId outgoingEdge = pairOfEdges.second;

		VertexId vertexStartRepeat = gp->g.EdgeEnd(incomingEdge); 
		VertexId vertexEndRepeat = gp->g.EdgeStart(outgoingEdge); 
	
		std::set<EdgeId> visited;
		dfs(vertexStartRepeat, vertexEndRepeat, component, visited, path);

		path.PushFront(incomingEdge);
		path.PushBack(outgoingEdge);

	
		// here we filter f.e. paths-loops which can be collected when performing dfs because there components go through 
		// vertices which are assigned with edges-loops so that they can be added in the dfs queue 
		if ( path.Size() < 3 ){
			return;
		}
	

	}

	void visit( const EdgeId& edge, std::set<EdgeId>& visited_edges, std::set<VertexId>& grey_vertices, 
		std::vector<EdgeId>& path, const std::vector<EdgeId>& components, 
		const std::vector<EdgeId>& singles, std::vector<EdgeId>& incoming_edges, std::vector<EdgeId>& outgoing_edges, bool& if_loop ) {

		VertexId edgeStartVertex = gp->g.EdgeStart(edge);
		VertexId edgeEndVertex = gp->g.EdgeEnd(edge);
		
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
		
			VertexId singleStartVertex = gp->g.EdgeStart(*single_edge);
			VertexId singleEndVertex = gp->g.EdgeEnd(*single_edge);

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

			VertexId componentEdgeStartVertex = gp->g.EdgeStart(*component_edge);
			VertexId componentEdgeEndVertex = gp->g.EdgeEnd(*component_edge);
			
			if ( componentEdgeStartVertex == edgeEndVertex || componentEdgeEndVertex == edgeStartVertex ) {
				
				//std::cout << "visiting edge " << gp->g.int_id(*component_edge) << std::endl;
				visit(*component_edge, visited_edges, grey_vertices, path, components, singles, incoming_edges, outgoing_edges, if_loop);
			}


		}

		grey_vertices.insert(edgeStartVertex);

	}


	void GetOtherEdges(path_extend::PathContainer& paths, const std::set<EdgeId>& used_edges){
	// adds edges from the rest of the graph (which was n)
		std::set<EdgeId> included;
		for (auto iter = gp->g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (used_edges.find(*iter) == used_edges.end() && included.find(*iter) == included.end()){
				paths.AddPair(new path_extend::BidirectionalPath(gp->g, *iter), new path_extend::BidirectionalPath(gp->g, gp->g.conjugate(*iter)));
				included.insert(*iter);
				included.insert(gp->g.conjugate(*iter));
			}
		}

	}


	void WriteResolved( const std::set<EdgeId>& used_edges, path_extend::PathContainer& resolved_paths, const std::string &file_name  ){

		path_extend::ContigWriter cw( gp->g);
		//cw.writePaths( resolvedPaths, fileName );
		//PathContainer paths;
		GetOtherEdges( resolved_paths, used_edges );
		cw.writePaths( resolved_paths, file_name );
	}
};
}

#endif
