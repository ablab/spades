/*
 * config_struct_hammer.cpp
 *
 *  Created on: Oct 15, 2011
 *      Author: snikolenko
 */

#include "config_struct_hammer.hpp"

void load(boost::property_tree::ptree const& pt, hammer_config& cfg)
{
	using config_common::load;
	// input options:
	load(pt, "working_dir", cfg.working_dir);
	load(pt, "reads", cfg.reads);
	load(pt, "tau", cfg.tau);
	load(pt, "num_threads", cfg.num_threads);
	load(pt, "num_iterations", cfg.num_iterations);
	load(pt, "quality_offset", cfg.quality_offset);

	load(pt, "paired_reads", cfg.paired_reads);
	if ( cfg.paired_reads ) {
		load(pt, "reads_left", cfg.reads_left);
		load(pt, "reads_right", cfg.reads_right);
	}

	load(pt, "read_blob_and_kmers", cfg.read_blob_and_kmers);
	load(pt, "write_blob_and_kmers", cfg.write_blob_and_kmers);
	load(pt, "exit_after_writing_blob_and_kmers", cfg.exit_after_writing_blob_and_kmers);
	if ( cfg.read_blob_and_kmers || cfg.write_blob_and_kmers ) {
		load(pt, "blob", cfg.blob);
		load(pt, "kmers", cfg.kmers);
	}

	load(pt, "read_kmers_after_clustering", cfg.read_kmers_after_clustering);
	load(pt, "write_kmers_after_clustering", cfg.write_kmers_after_clustering);
	load(pt, "kmers_after_clustering", cfg.kmers_after_clustering);

	load(pt, "error_rate", cfg.error_rate);
	load(pt, "blocksize_quadratic_threshold", cfg.blocksize_quadratic_threshold);
	load(pt, "good_cluster_threshold", cfg.good_cluster_threshold);
	load(pt, "blob_margin", cfg.blob_margin);
	load(pt, "trim_quality", cfg.trim_quality);

	load(pt, "trim_left_right", cfg.trim_left_right);
	load(pt, "use_iterative_reconstruction", cfg.use_iterative_reconstruction);
	load(pt, "reconstruction_in_full_iterations", cfg.reconstruction_in_full_iterations);
	load(pt, "iterative_reconstruction_threshold", cfg.iterative_reconstruction_threshold);
	load(pt, "max_reconstruction_iterations", cfg.max_reconstruction_iterations);
	load(pt, "write_each_iteration_kmers", cfg.write_each_iteration_kmers);
	load(pt, "regular_threshold_for_correction", cfg.regular_threshold_for_correction);
	load(pt, "discard_only_singletons", cfg.discard_only_singletons);
	load(pt, "special_nonsingleton_threshold", cfg.special_nonsingleton_threshold);
	load(pt, "use_true_likelihood", cfg.use_true_likelihood);

	load(pt, "conserve_memory", cfg.conserve_memory);
	load(pt, "num_of_tmp_files", cfg.num_of_tmp_files);
	load(pt, "skip_to_clustering", cfg.skip_to_clustering);
	load(pt, "skip_to_subvectors", cfg.skip_to_subvectors);
	load(pt, "skip_sorting_subvectors", cfg.skip_sorting_subvectors);
	load(pt, "unload_blob_before_merge", cfg.unload_blob_before_merge);

	load(pt, "likelihood_e_step", cfg.likelihood_e_step);
	load(pt, "subtract_simplex_volume", cfg.subtract_simplex_volume);

	load(pt, "debug_output_clustering", cfg.debug_output_clustering);
	load(pt, "debug_output_likelihood", cfg.debug_output_likelihood);
}
