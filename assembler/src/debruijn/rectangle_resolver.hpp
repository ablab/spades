/*
 * rectangle_resolver.hpp
 *
 *  Created on: Feb 7, 2012
 *      Author: Nikolay Vyahhi <vyahhi@gmail.com>
 */

#ifndef RECTANGLE_RESOLVER_HPP_
#define RECTANGLE_RESOLVER_HPP_

#include <string>
#include <log4cxx/logger.h>
#include "graph_pack.hpp"
#include "debruijn_stats.hpp"
#include "graphio.hpp"
#include <Python.h>

namespace debruijn_graph {

class RectangleResolver {
public:
	static void resolve(conj_graph_pack& gp) {
		INFO("Rectangle resolving started");
		const std::string output_dir = cfg::get().output_dir + "saves/";

		//conj_graph_pack resolved_gp(gp.genome);
		conj_graph_pack* resolved_gp = &gp;

		// Prepare input
		OutputContigs(resolved_gp->g, output_dir + "rectangle_before_contigs.fasta");
		PrintGraphPack(output_dir + "rectangle_before", *resolved_gp);

		// Init some parameters
		const std::string rr_filename = "src/debruijn/rr2.py";
		const std::string grp_filename = output_dir + "rectangle_before.grp";
		const std::string sqn_filename = output_dir + "rectangle_before.sqn";
		const std::string prd_filename = output_dir + "rectangle_before_et.prd";
		const std::string out_filename = output_dir + "rectangle_after";
		const std::string d = ToString(cfg::get().ds.IS.get() - cfg::get().ds.RL);

		// Ugly argc / argv
		int argc = 6;
		char *argv[6];
		for (int i = 0; i < 6; ++i) {
			argv[i] = new char[1000];
		}
		strcpy(argv[0], rr_filename.c_str());
		strcpy(argv[1], grp_filename.c_str());
		strcpy(argv[2], sqn_filename.c_str());
		strcpy(argv[3], prd_filename.c_str());
		strcpy(argv[4], d.c_str());
		strcpy(argv[5], out_filename.c_str());
		int updatepath = 1;

		// Run RR in Python
		INFO("Running rr2 in Python");
		Py_Initialize();
		PySys_SetArgvEx(argc, argv, updatepath);
		FILE* fp = fopen(rr_filename.c_str(), "r");
		PyRun_SimpleFile(fp, rr_filename.c_str());
		Py_Finalize();
		INFO("Saved to " + out_filename + ".grp .sqn _et.prd .fasta .log");

		// Free memory (somehow without this archaic things?)
		for (int i = 0; i < 6; ++i) {
			delete argv[i];
		}

		// Load RR output back from 'out_filename' saves
		// not implemented, TODO


		INFO("Rectangle resolving finished");
	}
};

}

#endif /* RECTANGLE_RESOLVER_HPP_ */
