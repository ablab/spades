#ifndef JOBWRAPERS_H_
#define JOBWRAPERS_H_

#include "running_modes.hpp"
#include "QcException.hpp"
#include "aho_corasick.hpp"
#include "output.hpp"
#include "config_struct_cclean.hpp"
#include "io/read_processor.hpp"
#include "ssw/ssw_cpp.h"

class AlignmentJobWrapper {
public:
	AlignmentJobWrapper(const Database * data, std::ostream& output, std::ostream& bed)
		:data(data), output(output), bed(bed){};

	bool operator()(const Read &r);

private:
	const Database * data;
	std::ostream& output;
	std::ostream& bed;
	const int mismatch_threshold = cclean_cfg::get().mismatch_threshold;
	const double aligned_part_fraction = cclean_cfg::get().aligned_part_fraction;
};

class ExactMatchJobWrapper {
public:
	ExactMatchJobWrapper(const Database * data, std::ostream& output, std::ostream& bed, const AhoCorasick &a)
		:data(data), output(output), bed(bed), ahoCorasick(a){};

	bool operator()(const Read &r);

private:
	const Database * data;
	std::ostream& output;
	std::ostream& bed;
	AhoCorasick ahoCorasick;
};

class ExactAndAlignJobWrapper {
public:
	ExactAndAlignJobWrapper(const Database * data, std::ostream& output, std::ostream& bed, const AhoCorasick &a)
		:data(data), output(output), bed(bed), ahoCorasick(a){};

	bool operator()(const Read &r);

private:
	const Database * data;
	std::ostream& output;
	std::ostream& bed;
	AhoCorasick ahoCorasick;
	const int mismatch_threshold = cclean_cfg::get().mismatch_threshold;
	const double aligned_part_fraction = cclean_cfg::get().aligned_part_fraction;
};

#endif /* JOBWRAPERS_H_ */
