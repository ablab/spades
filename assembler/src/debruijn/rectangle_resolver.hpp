/*
 * rectangle_resolver.hpp
 *
 *  Created on: Feb 7, 2012
 *      Author: Nikolay Vyahhi <vyahhi@gmail.com>
 */

#ifndef RECTANGLE_RESOLVER_HPP_
#define RECTANGLE_RESOLVER_HPP_

#include <string>

namespace debruijn_graph {

class RectangleResolver {
public:
	static void resolve(conj_graph_pack& gp) {
		INFO("Rectangle resolving started");

		const std::string output_dir = cfg::get().output_dir + "rectangle_resolve/";

		//conj_graph_pack resolved_gp(gp.genome);
		conj_graph_pack* resolved_gp = &gp;

		OutputContigs(resolved_gp->g, cfg::get().output_dir + "rectangle_contigs.fasta");
		PrintGraphPack(cfg::get().output_dir + "/saves/rectangle", *resolved_gp);
		// TODO: output

		INFO("Rectangle resolving finished");
	}
};

}

#endif /* RECTANGLE_RESOLVER_HPP_ */
