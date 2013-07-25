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

namespace debruijn_graph{

template <class GraphPack>
class CoverageBasedResolution {
	typedef double coverage_value;

	GraphPack *gp;
	std::vector< std::vector<EdgeId> > allPaths;

	//path with conjugate edges
	std::vector< std::vector<EdgeId> > resolvedPaths;

	const double threshold_one_list_;
	const double threshold_match_;
	const double threshold_global_;
	const double tandem_lower_threshold_;
	const double tandem_upper_threshold_;
	const double repeat_length_upper_threshold_;

	public:
	CoverageBasedResolution( GraphPack *gpack_arg, double threshold_one_list, double threshold_match,
				double threshold_global, double tandem_lower_threshold, double tandem_upper_threshold, double repeat_length_upper_threshold) :
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
	void joinPaths( std::vector<std::vector<EdgeId> >& paths, //std::vector< std::vector<EdgeId> >& paths,
			std::vector< std::vector<EdgeId> >& resolvedLoops,
			std::vector< std::vector<EdgeId> >& allPaths ) {

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

			INFO("path before");
			for ( auto e = path->begin(); e != path->end(); ++e ){
				std::cout << gp->g.int_id(*e) <<  " ";
			}
			std::cout << std::endl;
			bool updated = true;
			while (updated) {

				auto foundPath = startEdgeToPath.find( path->back() );
				if (foundPath != startEdgeToPath.end() ) {
					INFO("found path before");
					for ( auto e = foundPath->second.begin(); e != foundPath->second.end(); ++e ){
						std::cout << gp->g.int_id(*e) <<  " ";
					}
					std::cout << std::endl;
					path->insert(path->end(), foundPath->second.begin() + 1, foundPath->second.end());
				}

				else {
					updated = false;
				}

			}
			INFO("path after");
			for ( auto e = path->begin(); e != path->end(); ++e ){
				std::cout << gp->g.int_id(*e) <<  " ";
			}
			std::cout << std::endl;
			allPaths.push_back(*path);
		}
		INFO("out of path joining");

	}


	public :
	template <class DetailedCoverage, class EdgeQualityLabeler>
	void resolve_repeats_by_coverage( DetailedCoverage& coverage,
					size_t insert_size,
					EdgeLabelHandler<typename GraphPack::graph_t>& labels_after,
					const EdgeQualityLabeler& quality_labeler,
					PairedInfoIndexT<typename GraphPack::graph_t> & clustered_index,
					std::vector< PathInfo<typename GraphPack::graph_t> >& filteredPaths
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
		INFO("Getting components...");
		getComponents(  components, singles, labels_after, quality_labeler, clustered_index, filter.prohibitedEdges );
		//getComponentsWithReference(  components, singles, quality_labeler, prohibitedEdges );
		//getComponents( gp, components, singles, quality_labeler, unresolvedLoops );

		INFO("Traversing graph...");
		traverseComponents( components, singles, coverage, insert_size, resolvedPaths );

		std::set<EdgeId> usedEdges;

		std::cout << "Paths before joining: " << std::endl;
		for ( auto p = resolvedPaths.begin(); p != resolvedPaths.end(); ++p) {
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

		joinPaths(resolvedPaths, filter.resolvedLoops, allPaths);
		std::cout << "after joining: " << allPaths.size() << std::endl;
		filterConjugate( usedEdges, allPaths, filteredPaths);
		std::cout << "before filtering size " << allPaths.size() << " filtered size: " << filteredPaths.size() << std::endl;

		//FILE* file = fopen("/home/ksenia/path_resolved.log", "w");
		for ( auto p = allPaths.begin(); p != allPaths.end(); ++p) {
			//fprintf(file, "resolved path \n");
			for ( auto iter = p->begin(); iter != p->end(); ++iter ) {
				std::cout << gp->g.int_id(*iter) << " ";
				//fprintf(file, "%d ", gp->g.int_id(*iter));
				//fprintf(file, " ");
			}
			std::cout << std::endl;
			//fprintf(file,"\n");
		}

		std::cout << "-------------------------" << std::endl;
		for ( auto p = filteredPaths.begin(); p != filteredPaths.end(); ++p) {
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
		std::string fileName = cfg::get().output_dir + "resolved_by_coverage.fasta";
		//INFO("Writing result");

		path_extend::PathContainer pathsToOutput;
		for ( auto p = filteredPaths.begin(); p != filteredPaths.end(); ++p) {
			path_extend::BidirectionalPath* bidirectional_path = new path_extend::BidirectionalPath( gp->g );
			path_extend::BidirectionalPath* conjugate_path = new path_extend::BidirectionalPath( gp->g );
			auto tmpPath = p->getPath();
			for (auto it = tmpPath.begin(); it != tmpPath.end(); ++it ){

					bidirectional_path->PushBack(*it);
					EdgeId cedge = gp->g.conjugate(*it);
					conjugate_path->PushFront( cedge );
			}

			pathsToOutput.AddPair( bidirectional_path, conjugate_path );
		}

		WriteResolved( usedEdges, pathsToOutput, fileName);
	}

	private:

	void filterConjugate( std::set<EdgeId>& usedEdges,
				const std::vector< std::vector<EdgeId> > & paths,
				std::vector< PathInfo<typename GraphPack::graph_t> >& filteredPaths) {


		INFO("filtering conjugate edges");
		for ( auto path = paths.begin(); path != paths.end(); ++path) {

			bool ifInsert = true;
			for (auto e = path->begin(); e != path->end(); ++e) {
				if ( gp->g.conjugate(*e) == *e ) continue;
				if ( usedEdges.find(gp->g.conjugate(*e)) != usedEdges.end() ) {
//TODO:: this is not true, if we have autoreverse edge.
					ifInsert = false;
					break;
				}

			}
			if (ifInsert) {
				INFO("inserting");
				if (! gp->edge_pos.IsConsistentWithGenome(*path)) {
					std::cout << "not consistent with genome: ";
					for (auto iter = path->begin(); iter != path->end(); ++iter) {
						auto positions = gp->edge_pos.GetEdgePositions(*iter);
						std::cout << gp->g.int_id(*iter) << " (";
						for (auto pos = positions.begin(); pos != positions.end(); ++pos) {
							std::cout << pos->start() << " - " << pos->end() << " ";
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

		INFO("inserting paths into set of the used edges");
		for ( auto path = filteredPaths.begin(); path != filteredPaths.end(); ++path ) {
			auto p = path->getPath();
			std::cout << "size of p: " << p.size() << std::endl;
			for (auto e = p.begin(); e != p.end(); ++e) {
				usedEdges.insert(gp->g.conjugate(*e));
			}
		}

		INFO("out of filtering");
	}

	template<class EdgeQualityLabeler>
	void getComponentsWithReference( std::vector<EdgeId>& components, std::vector<EdgeId>& singles,
					EdgeLabelHandler<typename GraphPack::graph_t>& labels_after,
					EdgeQualityLabeler& quality_labeler ){

		INFO("Finding Components With Paired Info");
		for (auto iter = gp->g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {

			if (quality_labeler.quality(*iter) > 1) {

				components.push_back(*iter);
			}

			else {

				singles.push_back(*iter);
			}
		}
		/*
		for ( auto single_edge = singles.begin(); single_edge != singles.end(); ){

			bool ifLoop = false;

			VertexId sFrom = gp->g.EdgeStart(*single_edge);
			VertexId sInto = gp->g.EdgeEnd(*single_edge);

			for ( auto component_edge = components.begin(); component_edge != components.end(); ++component_edge) {


				VertexId cFrom = gp->g.EdgeStart(*component_edge);
				VertexId cInto = gp->g.EdgeEnd(*component_edge);

				if ( cFrom == sInto && sFrom == cInto ){

					components.push_back(*single_edge);
					singles.erase(single_edge);
					ifLoop = true;
					break;
				}
			}

			if (!ifLoop){
			 	++single_edge;
			}

		} */

		INFO("Components Ready");
/*
		std::cout << "SINGLES: ";
		for (auto sit = singles.begin(); sit != singles.end(); ++sit){
			std::cout << gp->g.int_id(*sit) << ", ";
		}
		std::cout << std::endl;

		std::cout << "COMPONENTS: ";
		for (auto cit = components.begin(); cit != components.end(); ++cit){
			std::cout << gp->g.int_id(*cit) << ", ";
		}
		std::cout << std::endl; */

	}

	bool isBulge( const EdgeId& edge ){

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
	bool checkIfComponentByPairedInfo( EdgeId edge, PairedInfoIndexT<Graph>& clustered_index, std::set<EdgeId>& prohibitedEdges ) {

//		auto improver = PairInfoImprover<Graph>(gp->g, clustered_index);
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
//				if (!improver.IsConsistent(edge, e1, e2, p1, p2)) {
//					std::cout << "Inconsistent for " << gp->g.int_id(edge) << ": " << gp->g.int_id(e1) << " " << gp->g.int_id(e2) << std::endl;
//					return true;
//				}
			}
		}

		return false;
	}

	template<class EdgeQualityLabeler>
	void getComponents( std::vector<EdgeId>& components, std::vector<EdgeId>& singles,
				EdgeLabelHandler<typename GraphPack::graph_t>& labels_after,
				const EdgeQualityLabeler& quality_labeler,
				PairedInfoIndexT<Graph>& clustered_index,
				std::set<EdgeId>& prohibitedEdges ){

		typedef int times;
		std::map<VertexId, times> out_degree;
		std::map<VertexId, times> in_degree;

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

		//TODO LengthCutoff move to config
		//double LengthCutoff = 0.0;
		/*int numEdges = 0;
		for (auto e_iter = gp->g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
			LengthCutoff += gp->g.length(*e_iter);
			numEdges += 1;
		}

		std::cout << "Length cutoff: " << 0.25 * LengthCutoff / numEdges << std::endl;
		*/
		//LengthCutoff = 300;//
		for (auto e_iter = gp->g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){

			if ( prohibitedEdges.find(*e_iter) != prohibitedEdges.end() ){
				continue;
			}

			VertexId from = gp->g.EdgeStart(*e_iter);
			VertexId into = gp->g.EdgeEnd(*e_iter);

			//std::cout << e_iter->int_id() << std::endl;
			if ( gp->g.length(*e_iter) >= cfg::get().rr.max_repeat_length ) {
				singles.push_back(*e_iter);
			}
			else if (! ( ((in_degree.find(from) == in_degree.end()) || out_degree[from] > 1) && ((out_degree.find(into) == out_degree.end()) || in_degree[into] > 1) )
			//else if (! ( ((in_degree.find(from) == in_degree.end()) || in_degree[from] > 1) || ((out_degree.find(into) == out_degree.end()) || out_degree[into] > 1) )
				) {

				components.push_back(*e_iter);
			}
			/*else if ( in_degree.find(from) != in_degree.end() && out_degree [from] == 2 && (out_degree.find(into) == out_degree.end()) || in_degree[into] > 1 ) {

					auto outgoingEdges = gp->g.OutgoingEdges(from);
					for ( auto e = outgoingEdges.begin(); e != outgoingEdges.end(); ++e ) {
						if ( prohibitedEdges.find(*e) != prohibitedEdges.end() )
							components.push_back(*e_iter);
							toComps = true;
					}

			}*/
			/*else if ( checkIfComponentByPairedInfo(*e_iter, clustered_index, prohibitedEdges ) ) {
				components.push_back(*e_iter);
				std::cout << "Component Edge Detected By Paired Info: " << gp->g.int_id(*e_iter) << std::endl;
			}*/

			// continue check if a short _terminal_ vertex
			else if( gp->g.length(*e_iter) < repeat_length_upper_threshold_ && (in_degree.find(from) != in_degree.end()) && (out_degree.find(into) != out_degree.end()) ){
			//	if (gp->g.int_id(*e_iter) ==10021207) INFO(2);

				components.push_back(*e_iter);
			}
			/*else if(gp->g.length(*e_iter) > 30000 ){
				singles.push_back(*e_iter);
			}*/
			else if ( labels_after.edge_inclusions.find(*e_iter) != labels_after.edge_inclusions.end() && labels_after.edge_inclusions[*e_iter].size()> 1)
			{
			//	if (gp->g.int_id(*e_iter) == 10021207) INFO(3);
				components.push_back(*e_iter);
			}
			else{
				singles.push_back(*e_iter);
			}

		}

		/* for (auto e_iter = gp->g.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter){
		 	EdgeId edgeid= *e_iter;
			std::cout << edgeid.int_id() << std::endl;
		 }

		INFO("AFTER GETTING COMPONENTS"); */

		/*for (auto degree_it = out_degree.begin(); degree_it != out_degree.end(); ++degree_it){
			std::cout << degree_it->first.int_id() << " " << degree_it->second << ", ";
		}
		std::cout << std::endl; */

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


	/*template <class GraphPack>
	bool containsSmallLoop( GraphPack& gp, const std::vector<EdgeId>& path){

		for ( auto it = path.begin(); it != path.end(); ++it ) {

			if (gp->g.EdgeStart(*it) == gp.g.EdgeEnd(*it)) {
				return true;
			}
		}
		return false;
	}*/

	bool containsSelfLoop( const std::vector<EdgeId>& path){

		for ( auto it = path.begin(); it != path.end(); ++it ) {

			if (gp->g.EdgeStart(*it) == gp->g.EdgeEnd(*it)) {
				return true;
			}
		}
		return false;
	}

	bool containsOnlyShortEdges( const std::vector<EdgeId>& path){

		for ( auto it = path.begin(); it != path.end(); ++it ) {

			if (gp->g.length(*it) >= repeat_length_upper_threshold_)
				return false;

		}
		return true;
	}


	void bfs ( const EdgeId& edge,  std::set<EdgeId>& visited_edges, const std::vector<EdgeId>& component, int& curLen, int& maxPathLen) {

		visited_edges.insert(edge);
		auto incomingEdges = gp->g.IncomingEdges(gp->g.EdgeStart(edge));

		for ( auto e = incomingEdges.begin(); e != incomingEdges.end(); ++e) {

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

	int getLongestPathLength( const std::vector<EdgeId>& component ){
	// gets a repetitive component and calculates the length of the longest path in it

		std::set<EdgeId> visited_edges;
		std::vector<std::vector<EdgeId>> paths;

		int maxPathLen = 0;
		for ( auto edge = component.begin(); edge != component.end(); ++edge ){

			if (visited_edges.find(*edge) != visited_edges.end()) continue;
			int curLen = gp->g.length(*edge);
			visited_edges.insert(*edge);
			bfs(*edge, visited_edges, component, curLen, maxPathLen);

		}

		return maxPathLen;

	}

	//todo WTF?!!!
	template <class DetailedCoverage>
	void traverseComponents( const std::vector<EdgeId>& components,
				const std::vector<EdgeId>& singles, DetailedCoverage& coverage,
				//std::set<EdgeId>& usedEdges,
				size_t insert_size,
				std::vector< std::vector<EdgeId>> & resolvedPaths) {

		std::set<EdgeId> visited_edges;
		INFO("Traversing components");
		FILE* file = fopen("/home/ksenia/path_resolved.log", "w");

		int allLength = 0;
		for (auto it = components.begin(); it != components.end(); ++it ){
			allLength += gp->g.length(*it);
		}

		int numberOfComponents = 0;
		int numberOfComponentWithDifferentInOutDegree = 0;

		int filteredByThresholds(0), resolvedPathsNum(0);
		for ( auto edge = components.begin(); edge != components.end(); ++edge ) {

			if ( visited_edges.find(*edge) != visited_edges.end() ){
				continue;
			}


			std::set<EdgeId> incomingEdges, outgoingEdges;
			std::vector<EdgeId> path;

			visit(*edge, visited_edges, path, components, singles, incomingEdges, outgoingEdges);


			int longestPathLen = getLongestPathLength(path);

			if ( containsSelfLoop(path) ) {
				continue;
			}

			if ( containsOnlyShortEdges(path) ) {
				INFO("contains only short edges");
				continue;
			}

			if ( insert_size < (size_t)longestPathLen ) numberOfComponents += 1;

			if (incomingEdges.size() != outgoingEdges.size() && insert_size < (size_t)longestPathLen )
				{

				numberOfComponentWithDifferentInOutDegree += 1;
				std::cout << "component with different in and out degree: " << longestPathLen << std::endl;
				for ( auto iter = path.begin(); iter != path.end(); ++iter ) {
					std::cout << gp->g.int_id(*iter)  << " ";
				}
				std::cout << std::endl;
				std::cout << "incoming edges: ";
				for ( auto iter = incomingEdges.begin(); iter != incomingEdges.end(); ++iter ) {
					std::cout << gp->g.int_id(*iter)  << " ";
				}
				std::cout << std::endl;
				std::cout << "outgoing edges: ";
				for ( auto iter = outgoingEdges.begin(); iter != outgoingEdges.end(); ++iter ) {
					std::cout << gp->g.int_id(*iter)  << " ";
				}
				std::cout << std::endl;
				continue;

			}
			if ( incomingEdges.size() == 0 || outgoingEdges.size() == 0) {
				//std::cout << incomingEdges.size() << " " << outgoingEdges.size() << std::endl;
				continue;
			}

			/*std::cout << "component: ";
			for ( auto iter = path.begin(); iter != path.end(); ++iter ) {
				std::cout << gp->g.int_id(*iter)  << " ";
			}
			std::cout << std::endl;
			*/

			//VERIFY_MSG(incomingEdges.size() != 0 and outgoingEdges.size() != 0, "vector of incoming and outgoing edges is zero length");
			std::vector<EdgeId> v(min(incomingEdges.size(), outgoingEdges.size()));
			std::vector<EdgeId>::iterator it = std::set_intersection(incomingEdges.begin(), incomingEdges.end(), outgoingEdges.begin(), outgoingEdges.end(), v.begin());

		 	v.resize(it - v.begin());

			/*if (v.size() > 0){

				std::cout << "LOOP:" << std::endl;
				for (auto it = v.begin(); it != v.end(); ++it){
					std::cout << gp->g.int_id(*it) << std::endl;
				}
				continue;
			} */
			fprintf(file, "resolved path \n");
			for ( auto iter = path.begin(); iter != path.end(); ++iter ) {
				fprintf(file, "%lu ", gp->g.int_id(*iter));
				fprintf(file, " ");
			}
			fprintf(file, "\nincoming edges: ");
			for ( auto e = incomingEdges.begin(); e != incomingEdges.end(); ++e) {

				fprintf(file, "%lu ", gp->g.int_id(*e));
			}
			fprintf(file,"\n");
			fprintf(file,"outgoing edges: ");
			for ( auto e = outgoingEdges.begin(); e != outgoingEdges.end(); ++e) {

				fprintf(file, "%lu ", gp->g.int_id(*e));
			}
			fprintf(file,"\n");

			fprintf(file,"\n");



			std::vector<std::pair<EdgeId, coverage_value>> incomingEdgesCoverage, outgoingEdgesCoverage;

			//std::cout << "still on" << std::endl;
			for ( auto inEdge = incomingEdges.begin(); inEdge != incomingEdges.end(); ++inEdge) {
				incomingEdgesCoverage.push_back(std::make_pair(*inEdge,coverage.GetOutCov(*inEdge)));
			}

			for ( auto outEdge = outgoingEdges.begin(); outEdge != outgoingEdges.end(); ++outEdge) {
				outgoingEdgesCoverage.push_back(std::make_pair(*outEdge,coverage.GetInCov(*outEdge)));
			}

			sort(incomingEdgesCoverage.begin(), incomingEdgesCoverage.end(), CompareSecond<EdgeId, coverage_value>());
			sort(outgoingEdgesCoverage.begin(), outgoingEdgesCoverage.end(), CompareSecond<EdgeId, coverage_value>());

		/*	INFO("incoming edges");
			for ( auto e = incomingEdgesCoverage.begin(); e != incomingEdgesCoverage.end(); ++e) {

				std::cout << e->first << " " << e->second << std::endl;
			}
			INFO("outgoing edges");
			for ( auto e = outgoingEdgesCoverage.begin(); e != outgoingEdgesCoverage.end(); ++e) {

				std::cout << e->first << " " << e->second << std::endl;
			}
		*/
			std::vector<std::pair<EdgeId,EdgeId>> pairsOfEdges;


			findClosest(incomingEdgesCoverage, outgoingEdgesCoverage, pairsOfEdges);
			//std::cout << "before repeat resolution " << incomingEdgesCoverage.size() << " " << outgoingEdgesCoverage.size() << " " << incomingEdges.size() << " " << path.size() << std::endl;
			if ( insert_size < (size_t)longestPathLen )
				if (pairsOfEdges.size() == 0)
					filteredByThresholds += 1;

			for ( auto edgePair = pairsOfEdges.begin(); edgePair != pairsOfEdges.end(); ++edgePair ){
				path_extend::BidirectionalPath* resolved_path = new path_extend::BidirectionalPath( gp->g );
				resolveRepeat( *edgePair, path, *resolved_path );
				if (resolved_path->Size() < 3 || resolved_path->Front() == resolved_path->Back() || coverage.GetOutCov(edgePair->first) < threshold_global_
					|| coverage.GetInCov(edgePair->second) < threshold_global_) {
					continue;
				}
				//now put it into collection of paths
				/*path_extend::BidirectionalPath* conjugate_path = new path_extend::BidirectionalPath( gp->g );
				bool skip = false;
				std::vector<EdgeId> tempPath = resolved_path->ToVector();
				for (auto it = tempPath.begin(); it != tempPath.end(); ++it ){

					EdgeId cedge = gp->g.conjugate(*it);
					if ( usedEdges.find(cedge) != usedEdges.end() ){
						skip = true;
						break;
					}
					conjugate_path->PushFront( cedge );
				} */

				bool skip = false;
				if (!skip) {

					std::vector<EdgeId> tempPath = resolved_path->ToVector();

					if ( insert_size < (size_t)longestPathLen )
						resolvedPathsNum += 1;
					resolvedPaths.push_back( tempPath );
					/*resolvedPaths.AddPair( resolved_path, conjugate_path );
					auto tmpPath = resolved_path->ToVector();
					if ( matchReference(*resolved_path, gp) ) {
						std::cout << "matches: ";
					}
					else {
						std::cout << "does not match: ";
					}
					for ( auto iter = tmpPath.begin(); iter != tmpPath.end(); ++iter ) {
						std::cout << gp->g.int_id(*iter) << " (" << coverage.inCoverage[*iter] << " - "<< coverage.outCoverage[*iter] << ") ";
					}
					std::cout << std::endl;

					for (auto it = tempPath.begin(); it != tempPath.end(); ++it ){
						usedEdges.insert( *it );
					}*/
				}
			}

		}

		std::cout << "Number of components: " << numberOfComponents << std::endl;
		std::cout << "Number of components with different in and out degree: " << numberOfComponentWithDifferentInOutDegree << std::endl;
		std::cout << filteredByThresholds << " " << resolvedPathsNum << std::endl;
	/*
		for ( auto iter = usedEdges.begin(); iter != usedEdges.end(); ++iter ){
			EdgeId cedge = gp->g.conjugate(*iter);
			usedEdges.insert(cedge);
		}
	*/
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

	void resolveRepeat( const std::pair<EdgeId,EdgeId>& pairOfEdges,
			const std::vector<EdgeId>& component, path_extend::BidirectionalPath& path ) {

		//std::cout << "resolveing repeat" << std::endl;
		EdgeId incomingEdge = pairOfEdges.first;
		EdgeId outgoingEdge = pairOfEdges.second;

		VertexId vertexStartRepeat = gp->g.EdgeEnd(incomingEdge);
		VertexId vertexEndRepeat = gp->g.EdgeStart(outgoingEdge);

		//std::cout << gp->g.int_id(incomingEdge) << " " << gp.g.int_id(outgoingEdge) << std::endl;
		std::set<EdgeId> visited;
		dfs(vertexStartRepeat, vertexEndRepeat, component, visited, path);

		path.PushFront(incomingEdge);
		path.PushBack(outgoingEdge);


		// here we filter f.e. paths-loops which can be collected when performing dfs because there components go through
		// vertices which are assigned with edges-loops so that they can be added in the dfs queue
		if ( path.Size() < 3 ){
			return;
		}

		/*std::cout << "PATH: ";
		for ( auto e = path.ToVector().begin(); e != path.ToVector().end(); ++e ){
			std::cout << gp->g.int_id(*e) << " ";
		}
		std::cout << std::endl;
		*/

	}

	void visit( const EdgeId& edge, std::set<EdgeId>& visited_edges,
		std::vector<EdgeId>& path, const std::vector<EdgeId>& components,
		const std::vector<EdgeId>& singles, std::set<EdgeId>& incoming_edges, std::set<EdgeId>& outgoing_edges ) {

		if (visited_edges.find(edge) != visited_edges.end()){
			return;
		}

		visited_edges.insert(edge);
		path.push_back(edge);
		VertexId edgeStartVertex = gp->g.EdgeStart(edge);
		VertexId edgeEndVertex = gp->g.EdgeEnd(edge);

		//std::cout << "visiting edge " << gp->g.int_id(edge) << " " << gp.g.int_id( edgeStartVertex ) << "," << gp.g.int_id( edgeEndVertex ) << std::endl;
		//VERIFY_MSG(singles.size() != 0, "size of singles in visit is zero!");
		if (singles.size() == 0) return;
		//std::cout << "vertices of single edges: ";
		for ( auto single_edge = singles.begin(); single_edge != singles.end(); ++single_edge ){

			VertexId singleStartVertex = gp->g.EdgeStart(*single_edge);
			VertexId singleEndVertex = gp->g.EdgeEnd(*single_edge);

			//std::cout << gp->g.int_id( *single_edge ) << " " << gp.g.int_id( singleStartVertex ) << " " << gp.g.int_id( singleEndVertex ) << std::endl;
			if ( singleStartVertex == edgeEndVertex ) {
				outgoing_edges.insert(*single_edge);
			}
			if ( singleEndVertex == edgeStartVertex ) {
				incoming_edges.insert(*single_edge);
			}
		}

		//std::cout << "SIZE OF OUTGOING AND INCOMING EDGES RESPECTIVELY " << outgoing_edges.size() << " " << incoming_edges.size() << std::endl;
		for ( auto component_edge = components.begin(); component_edge != components.end(); ++component_edge ){

			if ( gp->g.int_id(*component_edge) == gp->g.int_id(edge) ){
				continue;
			}
			//std::cout << "New component edge " << gp->g.int_id( *component_edge ) << std::endl;
			VertexId componentEdgeStartVertex = gp->g.EdgeStart(*component_edge);
			VertexId componentEdgeEndVertex = gp->g.EdgeEnd(*component_edge);

			if ( componentEdgeStartVertex == edgeEndVertex || componentEdgeEndVertex == edgeStartVertex ) {
				 visit(*component_edge, visited_edges, path, components, singles, incoming_edges, outgoing_edges);
			}

		}

		/*std::cout << "INCOMING ";
		for ( auto iterator = incoming_edges.begin(); iterator != incoming_edges.end(); ++iterator) {

			std::cout << gp->g.int_id( *iterator ) << " ";

		}

		std::cout << "OUTGOING ";
		for ( auto iterator = outgoing_edges.begin(); iterator != outgoing_edges.end(); ++iterator) {

			std::cout << gp->g.int_id( *iterator ) << " ";

		}

		std::cout << std::endl;*/

	}


	void getOtherEdges(path_extend::PathContainer& paths, const std::set<EdgeId>& usedEdges){
	// adds edges from the rest of the graph (which was n)
		std::set<EdgeId> included;
		for (auto iter = gp->g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
			if (usedEdges.find(*iter) == usedEdges.end() && included.find(*iter) == included.end()){
				paths.AddPair(new path_extend::BidirectionalPath(gp->g, *iter), new path_extend::BidirectionalPath(gp->g, gp->g.conjugate(*iter)));
				included.insert(*iter);
				included.insert(gp->g.conjugate(*iter));
			}
		}

	}


	void WriteResolved( const std::set<EdgeId>& usedEdges, path_extend::PathContainer& resolvedPaths, const std::string &fileName  ){

		path_extend::ContigWriter cw( gp->g);
		//cw.writePaths( resolvedPaths, fileName );
		//PathContainer paths;
		getOtherEdges( resolvedPaths, usedEdges );
		cw.writePaths( resolvedPaths, fileName );
	}
};
}

#endif
