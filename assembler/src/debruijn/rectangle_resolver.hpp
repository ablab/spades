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

		const std::string output_dir = cfg::get().output_dir + "rectangle_resolve/";

		//conj_graph_pack resolved_gp(gp.genome);
		conj_graph_pack* resolved_gp = &gp;

		// Prepare input
		OutputContigs(resolved_gp->g, output_dir + "before_rectangle_contigs.fasta");
		PrintGraphPack(output_dir + "/saves/before_rectangle", *resolved_gp);

		// Run RR in Python
		Py_Initialize();
		const std::string rr_filename = "src/debruijn/rr2.py";
		FILE* fp = fopen(rr_filename.c_str(), "r");
		PyRun_SimpleFile(fp, rr_filename.c_str());
		Py_Finalize();

		// Load RR output back

		INFO("Rectangle resolving finished");
	}
};

}

#endif /* RECTANGLE_RESOLVER_HPP_ */
