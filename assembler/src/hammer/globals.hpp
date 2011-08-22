#ifndef HAMMER_GLOBALS_HPP_
#define HAMMER_GLOBALS_HPP_

struct Globals {
	static int qvoffset;
	static double error_rate;
	static int blocksize_quadratic_threshold;
	static double good_cluster_threshold;
	static double blob_margin;
	static bool paired_reads;
	static int trim_quality;
};

#endif //  HAMMER_GLOBALS_HPP_

