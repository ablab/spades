/*
 * config_struct_hammer.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: snikolenko
 */

#include "config_struct_hammer.hpp"

void load(hammer_config& cfg, boost::property_tree::ptree const& pt)
{
	using config_common::load;
	// input options:
	load(cfg.working_dir   , pt, "working_dir"      );
	load(cfg.reads         , pt, "reads"            );
	load(cfg.tau           , pt, "tau"              );
	load(cfg.num_threads   , pt, "num_threads"      );
	load(cfg.num_iterations, pt, "num_iterations"   );
	load(cfg.quality_offset, pt, "quality_offset"   );

	load(cfg.paired_reads  , pt, "paired_reads"     );

    if ( cfg.paired_reads ) {
		load(cfg.reads_left , pt, "reads_left" );
		load(cfg.reads_right, pt, "reads_right");
	}

	load(cfg.read_blob_and_kmers              , pt, "read_blob_and_kmers"              );
	load(cfg.write_blob_and_kmers             , pt, "write_blob_and_kmers"             );
	load(cfg.exit_after_writing_blob_and_kmers, pt, "exit_after_writing_blob_and_kmers");


	if ( cfg.read_blob_and_kmers || cfg.write_blob_and_kmers ) {
		load(cfg.blob , pt, "blob" );
		load(cfg.kmers, pt, "kmers");
	}

	load(cfg.read_kmers_after_clustering ,  pt, "read_kmers_after_clustering" );
	load(cfg.write_kmers_after_clustering,  pt, "write_kmers_after_clustering");
	load(cfg.kmers_after_clustering      ,  pt, "kmers_after_clustering"      );

	load(cfg.error_rate                   , pt, "error_rate"                   );
	load(cfg.blocksize_quadratic_threshold, pt, "blocksize_quadratic_threshold");
	load(cfg.good_cluster_threshold       , pt, "good_cluster_threshold"       );
	load(cfg.blob_margin                  , pt, "blob_margin"                  );
	load(cfg.trim_quality                 , pt, "trim_quality"                 );

	load(cfg.trim_left_right                   , pt, "trim_left_right"                   );
	load(cfg.use_iterative_reconstruction      , pt, "use_iterative_reconstruction"      );
	load(cfg.reconstruction_in_full_iterations , pt, "reconstruction_in_full_iterations" );
	load(cfg.iterative_reconstruction_threshold, pt, "iterative_reconstruction_threshold");
	load(cfg.max_reconstruction_iterations     , pt, "max_reconstruction_iterations"     );
	load(cfg.write_each_iteration_kmers        , pt, "write_each_iteration_kmers"        );
	load(cfg.regular_threshold_for_correction  , pt, "regular_threshold_for_correction"  );
	load(cfg.discard_only_singletons           , pt, "discard_only_singletons"           );
	load(cfg.special_nonsingleton_threshold    , pt, "special_nonsingleton_threshold"    );
	load(cfg.use_true_likelihood               , pt, "use_true_likelihood"               );

	load(cfg.conserve_memory          ,pt, "conserve_memory"         );
	load(cfg.num_of_tmp_files         ,pt, "num_of_tmp_files"        );
	load(cfg.skip_to_clustering       ,pt, "skip_to_clustering"      );
	load(cfg.skip_to_subvectors       ,pt, "skip_to_subvectors"      );
	load(cfg.skip_sorting_subvectors  ,pt, "skip_sorting_subvectors" );
	load(cfg.skip_cluster_merging     ,pt, "skip_cluster_merging"    );
	load(cfg.skip_iterative           ,pt, "skip_iterative"          );
	load(cfg.hamming_class_buffer     ,pt, "hamming_class_buffer"    );
	load(cfg.unload_blob_before_merge ,pt, "unload_blob_before_merge");

	load(cfg.likelihood_e_step      , pt, "likelihood_e_step"      );
	load(cfg.subtract_simplex_volume, pt, "subtract_simplex_volume");
	load(cfg.change_n_to_random     , pt, "change_n_to_random"     );

	load(cfg.debug_output_clustering, pt, "debug_output_clustering");
	load(cfg.debug_output_likelihood, pt, "debug_output_likelihood");

	// TODO: remove old options when done

	load(cfg.general_do_everything_after_first_iteration, pt, "general_do_everything_after_first_iteration");
	load(cfg.general_reconstruct_only, pt, "general_reconstruct_only");
	load(cfg.general_change_n_to_a, pt, "general_change_n_to_a");
	load(cfg.general_hard_memory_limit, pt, "general_hard_memory_limit");
	load(cfg.general_max_nthreads, pt, "general_max_nthreads");
	load(cfg.general_tau, pt, "general_tau");
	load(cfg.general_max_iterations, pt, "general_max_iterations");
	load(cfg.general_blob_margin, pt, "general_blob_margin");

	load(cfg.count_do, pt, "count_do");
	load(cfg.count_numfiles, pt, "count_numfiles");
	load(cfg.count_merge_nthreads, pt, "count_merge_nthreads");

	load(cfg.sort_do, pt, "sort_do");

	load(cfg.subvectors_do, pt, "subvectors_do");
	load(cfg.subvectors_blocksize_quadratic_threshold, pt, "subvectors_blocksize_quadratic_threshold");

	load(cfg.hamming_do, pt, "hamming_do");
	load(cfg.hamming_write_bad_kmers, pt, "hamming_write_bad_kmers");
	load(cfg.hamming_write_solid_kmers, pt, "hamming_write_solid_kmers");

	load(cfg.bayes_do, pt, "bayes_do");
	load(cfg.bayes_nthreads, pt, "bayes_nthreads");
	load(cfg.bayes_quality_threshold, pt, "bayes_quality_threshold");
	load(cfg.bayes_singleton_threshold, pt, "bayes_singleton_threshold");
	load(cfg.bayes_nonsingleton_threshold, pt, "bayes_nonsingleton_threshold");
	load(cfg.bayes_discard_only_singletons, pt, "bayes_discard_only_singletons");
	load(cfg.bayes_debug_output, pt, "bayes_debug_output");
	load(cfg.bayes_use_hamming_dist, pt, "bayes_use_hamming_dist");

	load(cfg.expand_do, pt, "expand_do");
	load(cfg.expand_max_iterations, pt, "expand_max_iterations");
	load(cfg.expand_nthreads, pt, "expand_nthreads");
	load(cfg.expand_write_each_iteration, pt, "expand_write_each_iteration");
	load(cfg.expand_write_kmers_result, pt, "expand_write_kmers_result");

	load(cfg.correct_do, pt, "correct_do");
	load(cfg.correct_threshold, pt, "correct_threshold");
	load(cfg.correct_use_threshold, pt, "correct_use_threshold");

	load(cfg.input_numfiles, pt, "input_numfiles");
	if (cfg.input_numfiles > 0) load(cfg.input_file_0, pt, "input_file_0");
	if (cfg.input_numfiles > 1) load(cfg.input_file_1, pt, "input_file_1");
	if (cfg.input_numfiles > 2) load(cfg.input_file_2, pt, "input_file_2");
	if (cfg.input_numfiles > 3) load(cfg.input_file_3, pt, "input_file_3");
	if (cfg.input_numfiles > 4) load(cfg.input_file_4, pt, "input_file_4");
	load(cfg.input_paired, pt, "input_paired");
	load(cfg.input_solid_kmers, pt, "input_solid_kmers");
	load(cfg.input_working_dir, pt, "input_working_dir");
	load(cfg.input_trim_quality, pt, "input_trim_quality");
	load(cfg.input_qvoffset, pt, "input_qvoffset");
	load(cfg.input_read_solid_kmers, pt, "input_read_solid_kmers");
	load(cfg.input_solid_kmers, pt, "input_solid_kmers");
}
