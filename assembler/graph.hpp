/* de Bruijn graph
 * 
 * (it should be more optimized -- it's a draft!)
 */

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <map>
#include "read.hpp"

struct Node { // 4*1 + 4*2 = 12 bytes
	char out_edges[4];
	short out_count[4]; // may be char
};

class Graph { // de Bruijn
	public:
		Graph(const std::list<MateRead> &mr);
		std::pair<const Kmer,Node> get_node(const Kmer k); // &nodes_.find(k)
	private:
		std::map<Kmer,Node> nodes_; // should be some hash_map (default std::map is red-black tree).
};

#endif
