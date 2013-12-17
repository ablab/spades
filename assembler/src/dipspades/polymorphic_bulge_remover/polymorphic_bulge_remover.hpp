#pragma once

#include "simple_bulge_remover.hpp"
#include "complex_bulge_remover.hpp"
#include "equal_sequence_gluer.hpp"
#include "iterative_tails_gluing.hpp"
#include "../../debruijn/stats/debruijn_stats.hpp"

#include "omni/visualization/visualization.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/graph_component.hpp"

using namespace debruijn_graph;

namespace dipspades {

class PolymorphicBulgeRemoverHelper {
public:
	typedef ComplexBulgeGluer<RelatedBaseGlueDirectionDefiner, GluingVericesDefiner, BulgeSplitter> BaseBulgeGluer;
	static BaseBulgeGluer CreateBaseBulgeGluer(Graph &graph, double rel_len_threshold){
		return BaseBulgeGluer(graph, RelatedBaseGlueDirectionDefiner(graph),
				GluingVericesDefiner(graph, rel_len_threshold), BulgeSplitter(graph));
	}
};

class PolymorphicBulgeRemover {
	conj_graph_pack 		&graph_pack_;
	BaseHistogram<size_t> 	&bulge_len_hist_;
	Compressor<Graph> 		compressor_;

	typedef BulgeRemoverAlgorithm<DijkstraBulgePathsSearcher,
			PolymorphicBulgeRemoverHelper::BaseBulgeGluer> LightBulgeRemover;
	typedef BulgeRemoverAlgorithm<PathProcessorBulgeSearcher,
			PolymorphicBulgeRemoverHelper::BaseBulgeGluer> HardBulgeRemover;

	void RunSimpleBRCycle(){
		INFO("Simple polymorphic bulge remover runs");
		SimpleBulgeRemover spath_br(graph_pack_.g, bulge_len_hist_, dsp_cfg::get().pbr);
		size_t num_glued_bulges = 1;
		for(size_t num_iter = 1; num_glued_bulges > 0; num_iter++){
			num_glued_bulges = spath_br.Run();
			compressor_.CompressAllVertices();
			INFO(ToString(num_iter) + " iteration: " + ToString(num_glued_bulges) + " simple bulges were glued");
		}
		INFO("Simple polymorphic bulge remover ends");
	}

	template<class BulgeRemover>
	void BulgeRemoverCycle(string bulge_remover_name, size_t num_iters){
		INFO(bulge_remover_name + " starts");
		INFO("Maximal number of iterations: " << num_iters);
		BulgeRemover br(graph_pack_.g,
				PolymorphicBulgeRemoverHelper::CreateBaseBulgeGluer(graph_pack_.g,
						dsp_cfg::get().pbr.paired_vert_rel_threshold),
				bulge_len_hist_,
				dsp_cfg::get().pbr);
		size_t num_glued_bulges = 1;
		for(size_t i = 0; (i < num_iters) && (num_glued_bulges != 0); i++){
			num_glued_bulges = br.Run();
			compressor_.CompressAllVertices();
			INFO(ToString(i + 1) + " iteration: " + ToString(num_glued_bulges) + " complex bulges were glued");
		}
		INFO(bulge_remover_name + " ends");
	}

	void WriteComponents(string component_dir) {
		if(!dsp_cfg::get().rp.developer_mode)
			return;

		make_dir(dsp_cfg::get().io.output_dir + "components/");
	    omnigraph::DefaultLabeler<Graph> labeler(graph_pack_.g, graph_pack_.edge_pos);
	    make_dir(dsp_cfg::get().io.output_dir + "components/" + component_dir + "/");
        omnigraph::visualization::WriteComponents(graph_pack_.g,
        		dsp_cfg::get().io.output_dir + "components/" + component_dir + "/",
        		omnigraph::ReliableSplitter<Graph>(graph_pack_.g),
        		omnigraph::visualization::DefaultColorer(graph_pack_.g, Path<EdgeId>(), Path<EdgeId>()),
        		labeler);
	}

public:
	PolymorphicBulgeRemover(conj_graph_pack &graph_pack,
			BaseHistogram<size_t> &bulge_len_hist) :
		graph_pack_(graph_pack),
		bulge_len_hist_(bulge_len_hist),
		compressor_(graph_pack_.g, false)  { }

	void Run(){
		if(!dsp_cfg::get().pbr.enabled)
			return ;

		WriteComponents("before_pbr");

		graph_pack_.kmer_mapper.SetUnsafeMode(true);

		INFO("Polymorphic bulge remover starts");
		RunSimpleBRCycle();
		BulgeRemoverCycle<LightBulgeRemover>("LightBulgeRemover", dsp_cfg::get().pbr.num_iters_lbr);

//		INFO("Iterative tail gluing starts");
//		IterativeTailGluing(graph_pack_.g, dsp_cfg::get().pbr.rel_bulge_align).IterativeProcessTails();
//		compressor_.CompressAllVertices();
//		INFO("Iterative tail gluing ends");

		INFO("Index refilling");
		graph_pack_.index.Refill();

        INFO("Glueing equal kmers");
        EqualSequencesGluer<Graph>(graph_pack_.g, graph_pack_.index).GlueEqualKmers();
		INFO("Polymorphic ends remover ends");

		WriteComponents("after_pbr");

		graph_pack_.kmer_mapper.SetUnsafeMode(false);
	}
};

}
