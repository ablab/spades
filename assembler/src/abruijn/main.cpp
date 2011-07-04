#include "toyexamples.hpp"
#include "graphBuilder.hpp"

#include "logging.hpp"
DECL_PROJECT_LOGGER("a")

#include "tip_clipper.hpp"
#include "bulge_remover.hpp"
#include "erroneous_connection_remover.hpp"
#include "statistics.hpp"
#include "omni_tools.hpp"
#include "omnigraph.hpp"
#include "visualization_utils.hpp"
#include "libs/getopt_pp/getopt_pp_standalone.h"
#include <iostream>

using namespace std;
using namespace GetOpt;

namespace abruijn {

DECL_LOGGER("main")

typedef omnigraph::Omnigraph Graph;

class Launch {
public:
	int take_;
	bool output_single_;
	size_t cut_;
	int mode_;

	int tc_max_coverage_;
	double tc_max_relative_coverage_;
	int tc_max_tip_length_;

	double br_max_coverage_;
	double br_max_relative_coverage_;
	double br_max_delta_;
	double br_max_relative_delta_;
	size_t br_max_length_div_K_;

	double ecr_max_coverage_;
	int ecr_max_length_div_K_;

	Graph* g_;

private:
	GetOpt_pp& options_;
	bool help_;
	bool log_;
	stringstream usage_;
	string input_files_[2];
	string output_file_;

	string option(char short_opt, string long_opt, string help_msg) {
		return option<string>(short_opt, long_opt, help_msg, "");
	}

	template<typename T>
	T option(char short_opt, string long_opt, string help_msg, T default_value, bool show_default_value = false) {
		T t;
		options_ >> Option(short_opt, long_opt, t, default_value);
		if (help_msg != "") {
			usage_ << "--" << long_opt << " " << help_msg;
			if (show_default_value) {
				usage_ << " (default = " << default_value << ")";
			}
			usage_ << "\n";
		}
		if (log_) {
			INFO(long_opt << " = " << t);
		}
		return t;
	}

	bool optionPresent(char short_opt, string long_opt, string help_msg) {
		bool t;
		options_ >> OptionPresent(short_opt, long_opt, t);
		if (help_msg != "") {
			usage_ << "--" << long_opt << " (with no value) " << help_msg << "\n";
		}
		if (log_) {
			INFO(long_opt << " = " << (t ? "true" : "false"));
		}
		return t;
	}

public:
	Launch(GetOpt::GetOpt_pp& options) : options_(options) {
		log_ = false;
		usage_ << "\nSet option values like this: \"--key1 value1 --key2 value2 ...\".\n";
		usage_ << "\nSupported options are:\n";

		input_files_[0] = input_files_[1] = option('i', "input", "first data file (with path). Second file name is obtained by replacing last '1' with '2'");
		for (int i = input_files_[1].length() - 1; i >= 0; i--) {
			if (input_files_[1][i] == '1') {
				input_files_[1][i] = '2';
				break;
			}
		}

		help_ = optionPresent('h', "help", "") || (input_files_[0] == "");
		log_ = !help_;
		if (log_) {
			INFO("Hello, A Bruijn! Using parameters:");
			INFO("input = " << input_files_[0] << " and " << input_files_[1]);
		}

		output_file_ = option('o', "output", "output file (with path)");
		take_ = option('t', "take", "how many minimizers to take from each read", 1, true);
		output_single_ = optionPresent('s', "single", "whether to visualize vertices non-RC-paired");
		cut_ = option('c', "cut", "how many first reads from the input file to use", -1, true);

		mode_ = option('m', "mode", "", 0);
		usage_ << "--mode bit mask (default = 0):\n";
		usage_ << "       1 whether to find minimizers, not local minimizers\n";
		usage_ << "       2 whether to ensure at least 2 minimizers in each read\n";
		if (log_) {
			if (mode_ & 1) {
				INFO("mode => find minimizers");
			} else {
				INFO("mode => find local minimizers");
			}
			if (mode_ & 2) {
				INFO("mode => ensure at least two minimizers in each read");
			}
		}

		usage_ << "\nTip clipping parameters:\n";
		tc_max_coverage_ = option(' ', "tc-max-coverage", " ", 1000, true);
		tc_max_relative_coverage_ = option(' ', "tc-max-relative-coverage", " ", 2.0, true);
		tc_max_tip_length_ = option(' ', "tc-max-tip-length", " ", 50, true);

		usage_ << "\nBulge removing parameters:\n";
		br_max_coverage_ = option(' ', "br-max-coverage", " ", 1000, true);
		br_max_relative_coverage_ = option(' ', "br-max-relative-coverage", " ", 1.1, true);
		br_max_delta_ = option(' ', "br-max-delta", " ", 4.0, true);
		br_max_relative_delta_ = option(' ', "br-max-relative-delta", " ", 0.1, true);
		br_max_length_div_K_ = option(' ', "br-max-length-div-k", " ", 5, true);

		usage_ << "\nErroneous connection removing parameters:\n";
		ecr_max_coverage_ = option(' ', "erc-max-coverage", " ", 20, true);
		ecr_max_length_div_K_ = option(' ', "erc-max-length-div-k", " ", 5, true);

		if (help_) {
			std::cout << usage_.str();
			exit(0);
		}
	}

	void stats(string name) {
		omnigraph::StrGraphLabeler<Graph> labeler(*g_);
		gvis::WriteToDotFile(output_file_ + "_" + name, "earmarked", *g_, labeler);
		INFO("Statistics of " << name << ":");
		omnigraph::VertexEdgeStat<Graph> vertex_edge_stat(*g_);
		vertex_edge_stat.Count();
	}

	void run() {
		StrobeReader<2, Read, ireadstream> sr(input_files_);
		PairedReader<ireadstream> paired_stream(sr, 220);
		SimpleReaderWrapper<PairedReader<ireadstream> > srw(paired_stream);
		CuttingReader<SimpleReaderWrapper<PairedReader<ireadstream> > > cr(srw, cut_);

		abruijn::GraphBuilderMaster<CuttingReader<SimpleReaderWrapper<PairedReader<ireadstream>>>> gbm(cr, take_, mode_);
		g_ = gbm.build();
		stats("uncompressed");

		//  INFO("Spelling the reference genome");
		//  gbm.SpellGenomeThroughGraph(cut_ + 219);

		INFO("===== Compressing... =====");
		omnigraph::Compressor<Graph> compressor(*g_);
		compressor.CompressAllVertices();
		stats("compressed");

		INFO("===== Clipping tips... =====");
		omnigraph::TipComparator<Graph> comparator(*g_);
		omnigraph::TipClipper<Graph, omnigraph::TipComparator<Graph>> tip_clipper(*g_, comparator, tc_max_tip_length_, tc_max_coverage_, tc_max_relative_coverage_);
		tip_clipper.ClipTips();
		stats("tc");

		INFO("===== Removing bulges... =====");
		omnigraph::SimplePathCondition<Graph> simple_path_condition(*g_);
		omnigraph::BulgeRemover<Graph, omnigraph::SimplePathCondition<Graph>> bulge_remover(*g_, br_max_length_div_K_ * K, br_max_coverage_, br_max_relative_coverage_, br_max_delta_, br_max_relative_delta_, simple_path_condition);
		bulge_remover.RemoveBulges();
		stats("tc_br");

		INFO("===== Removing erroneous connections... =====");
		omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(ecr_max_length_div_K_ * K, ecr_max_coverage_);
		erroneous_edge_remover.RemoveEdges(*g_);
		stats("tc_br_ecr");

		INFO("===== Clipping tips #2... =====");
		tip_clipper.ClipTips();
		stats("tc_br_ecr_tc");

		INFO("===== Removing bulges #2... =====");
		bulge_remover.RemoveBulges();
		stats("tc_br_ecr_tc_br");

//		INFO("===== Outputting graph to " << output_file_ << " =====");
//		gvis::WriteToDotFile(output_file_, "earmarked", g_, labeler);

		//ABruijnGraphWithGraphVisualizer ( "ATGTGTGACTTTGTATCGTATTGCGGGCGGCGCGCTTATTGTATGCGTAAATTTGGGTCATATTGATCGTAAAATGCGTATGATGCACTGCA", 6, 3 );
	}
};

} // namespace abruijn


int main(int argc, char* argv[]) {
	GetOpt_pp options(argc, argv, Include_Environment);
	abruijn::Launch launch(options);
	launch.run();
	return 0;
}
