#ifndef LOOP_FILTER
#define LOOP_FILTER

#include <algorithm>
#include "graph_pack.hpp"
#include "omni/mf_ec_remover.hpp"

namespace debruijn_graph {
	
	template <class graph_pack, class DetailedCoverage>
	class LoopFilter {
		const graph_pack& gp;
		const DetailedCoverage& coverage_;
		std::vector<VertexId> order;
		set<VertexId> used_vertices_;
		std::vector< set<EdgeId> > simple_loops_;
		//contains all the vertices that go to complex loops
		set<EdgeId> prohibited_edges_;
		std::vector< std::vector<EdgeId> > resolved_loops_;
		const double ratio_lower_threshold_;
		const double ratio_upper_threshold_;
		const double repeat_length_upper_threshold_;

		public:
		LoopFilter ( const graph_pack& g, const DetailedCoverage& index, double ratio_lower_threshold, double ratio_upper_threshold, double repeat_length_upper_threshold ):
													gp(g),
													coverage_(index),
													ratio_lower_threshold_(ratio_lower_threshold),
													ratio_upper_threshold_(ratio_upper_threshold),
													repeat_length_upper_threshold_(repeat_length_upper_threshold) {
		}


		const vector<vector<EdgeId>>& resolved_loops() const {
			return resolved_loops_;
		}

		void get_loopy_components( const EdgeQuality<typename graph_pack::graph_t>& quality_labeler ) {
			INFO("Detecting loops...");
			for ( auto v = gp.g.begin(); v != gp.g.end(); ++v ) {
				if ( used_vertices_.find(*v) != used_vertices_.end() ) 
					continue;
				dfs_down(*v);
			}
			used_vertices_.clear();
			std::vector<std::vector<VertexId>> loops;
			for ( unsigned i = 0; i < order.size(); ++i ) {
				VertexId startVertex = order[order.size() - i - 1];
				if ( used_vertices_.find(startVertex) != used_vertices_.end() ) 
					continue;
				std::vector<VertexId> loop;
				dfs_up(startVertex, loop);
				loops.push_back(loop);
			}
			int L = 0;
			int loops_counter = 0;
			int simple_loops_counter = 0;
			for (auto loop = loops.begin(); loop != loops.end(); ++loop ){
				if (loop->size() == 1) {
					continue;
				}
				loops_counter++;
				set<EdgeId> path;
				std::vector<EdgeId> resolved_loop;
				EdgeId incomingEdge, outgoingEdge;
				bool resolved = false;
				bool ifSimple = IfSimpleLoop(*loop, path, incomingEdge, outgoingEdge);
				if (ifSimple) {
					simple_loops_counter++;
					if ( gp.g.int_id(incomingEdge) != 0 && gp.g.int_id(outgoingEdge) != 0 )
						resolved = ResolveSimpleLoop( resolved_loop, incomingEdge, outgoingEdge );
				}
				/*if (!resolved) {
					prohibited_edges_.insert(path.begin(), path.end());
				}*/
				if (resolved) {
					for (auto e = resolved_loop.begin(); e != resolved_loop.end(); ++e ) {
						L += gp.g.length(*e);
					}
					resolved_loops_.push_back(resolved_loop);
				}

				string str = "";
				if (!resolved) str = "not ";

				DEBUG(str << "resolved ");
				for (auto e = resolved_loop.begin(); e != resolved_loop.end(); ++e ) {
					DEBUG(gp.g.int_id(*e) << "q: (" << quality_labeler.quality(*e) << ") ");
				}
				DEBUG("\n");
			}
			DEBUG("Overall length is " << L << ", number of loops: " << loops_counter << " " 
					<< "number of resolved loops: " << resolved_loops_.size() << " " 
						<< "simple loops: " << simple_loops_counter << "\n");

		}	
	
		private:
		bool ResolveSimpleLoop( std::vector<EdgeId>& resolved_loop, 
					const EdgeId& incomingEdge, const EdgeId& outgoingEdge ) const  {

			EdgeId startEdge, endEdge;
			bool can_be_resolved = true;
			VertexId inVertex = gp.g.EdgeEnd( incomingEdge );
			auto inCov = coverage_.GetOutCov(incomingEdge);
			auto outCov = coverage_.GetInCov(outgoingEdge);
			auto cov = (inCov + outCov) / 2.0;	
			DEBUG("inCoverage: " << inCov << "; outCoverage " << outCov << "; cov: " << cov << "\n");
			auto loopEdge1 = *gp.g.OutgoingEdges(inVertex).begin();
			VertexId edge1End = gp.g.EdgeEnd(loopEdge1);
			EdgeId loopEdge2;
			auto loopEdge2it = gp.g.OutgoingEdges(edge1End).begin();
			if ( *loopEdge2it == outgoingEdge ) {
				++loopEdge2it;
			}
			loopEdge2 = *loopEdge2it;
			auto cov1 = gp.g.coverage(loopEdge1);
			auto cov2 = gp.g.coverage(loopEdge2);
			double intpart;
			auto ratio1 = modf( cov1 / cov, &intpart );
			auto ratio2 = modf( cov2 / cov, &intpart );
			if ( ( ratio1 > ratio_lower_threshold_ && ratio1 < ratio_upper_threshold_ ) || ( ratio2 > ratio_lower_threshold_ && ratio2 < ratio_upper_threshold_ ) ) {
				can_be_resolved = false;
			}
			auto time1 = floor( cov1 / cov + 0.5 );
			auto time2 = floor( cov2 / cov + 0.5 );
			if (time1 - time2 != 1 || time2 == 0) {
				can_be_resolved = false;
			}
			//if (can_be_resolved) {
				resolved_loop.push_back(incomingEdge);
				resolved_loop.push_back(loopEdge1);
				for ( int i = 0; i < time2; ++i ) {
					resolved_loop.push_back(loopEdge2);
					resolved_loop.push_back(loopEdge1);
				}
				resolved_loop.push_back(outgoingEdge);
				if ( resolved_loop.size() > 5 ) can_be_resolved = false;
			//}
			return can_be_resolved;
		}

		bool IfSimpleLoop( vector<VertexId>& loop, set<EdgeId>& path, EdgeId& incomingEdge, EdgeId& outgoingEdge ) const {
		//the idea: if a loop is simple there is a single possible way to come into any vertex in the component
			bool ifSimple = true;
			int nextVerticesOutOfLoop = 0;
			int prevVerticesOutOfLoop = 0;
			for ( auto v = loop.begin(); v != loop.end(); ++v ) {
				auto incomingEdges = gp.g.IncomingEdges(*v);
				int prevVerticesInLoop = 0;
				for ( auto e = incomingEdges.begin(); e != incomingEdges.end(); ++e ){
					auto startVertex = gp.g.EdgeStart(*e);
					// count the number of edges in the loop coming into this vertex	
					if ( find(loop.begin(),loop.end(),startVertex) != loop.end() ){
						prevVerticesInLoop+=1;
						if (prevVerticesInLoop > 1) {
							ifSimple = false;
						}
						path.insert(*e);
						//prohibited_edges_.insert(*e);
					}
					//count the number of edges coming into loop, i.e. check if the loop is a repeat
					else {
						prevVerticesOutOfLoop+=1;
						incomingEdge = *e;
						if (prevVerticesOutOfLoop > 1) {
							ifSimple = false;
						}
					}
				}
				auto outgoingEdges = gp.g.OutgoingEdges(*v);
				int nextVerticesInLoop = 0;
				for ( auto e = outgoingEdges.begin(); e != outgoingEdges.end(); ++e ){
					auto endVertex = gp.g.EdgeEnd(*e);
					// count the number of edges in the loop coming out of the vertex	
					if ( find(loop.begin(), loop.end(), endVertex) != loop.end() ){
						nextVerticesInLoop += 1;
						if (nextVerticesInLoop > 1){
							ifSimple = false;
						}
						path.insert(*e);
						//prohibited_edges_.insert(*e);
					}
					//count the number of edges coming into loop, i.e. check if the loop is a repeat
					else {
						nextVerticesOutOfLoop+=1;
						outgoingEdge = *e;
						if (nextVerticesOutOfLoop > 1) {
							ifSimple = false;
						}
					}

				}

			}
			return ifSimple;
		}

		void dfs_down( const VertexId& v ){
			used_vertices_.insert(v);
			for ( auto incidentEdge = gp.g.out_begin(v); incidentEdge != gp.g.out_end(v); ++incidentEdge ){
				if ( gp.g.length(*incidentEdge) > cfg::get().rr.max_repeat_length  ) continue;
				VertexId vOut = gp.g.EdgeEnd(*incidentEdge);
				if (used_vertices_.find(vOut) != used_vertices_.end())
					continue;
				dfs_down( vOut );
			}
			order.push_back(v);
		}

		void dfs_up( const VertexId& v, vector<VertexId>& loop ){
			used_vertices_.insert( v );
			loop.push_back( v );
			for ( auto incidentEdge = gp.g.in_begin(v); incidentEdge != gp.g.in_end(v); ++incidentEdge ) {
				if ( gp.g.length(*incidentEdge) > cfg::get().rr.max_repeat_length  ) continue;
				VertexId vIn = gp.g.EdgeStart(*incidentEdge);
				if (used_vertices_.find(vIn) != used_vertices_.end()) 
					continue;
				dfs_up( vIn, loop );
			}
		}
	};
}

#endif
