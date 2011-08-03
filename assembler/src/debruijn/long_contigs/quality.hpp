/*
 * quality.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef QUALITY_HPP_
#define QUALITY_HPP_

namespace long_contigs {

using namespace debruijn_graph;


// ====== Quality functions ======
//Find bidirectional path in given genome path
size_t FindInGenomePath(BidirectionalPath& myPath, Path<Graph::EdgeId>& genomePath) {
	if (myPath.size() > genomePath.size()) {
		INFO("Warning, unexpected path length")
		return -1;
	}

	for (size_t i = 0; i < genomePath.size() - myPath.size() + 1; ++i) {
		bool found = true;

		size_t j = i;
		for(auto iter = myPath.begin(); iter != myPath.end(); ++iter) {
			if (*iter != genomePath[j]) {
				found = false;
				break;
			}
			++j;
		}

		if (found) {
			return i;
		}
	}
	return -1;
}

//Find inexact match to genome path
size_t FindInGenomeInexact(Graph& g, BidirectionalPath& myPath, Path<Graph::EdgeId>& genomePath, int& startPos, size_t& maxLengthMached) {
	if (myPath.size() > genomePath.size()) {
		INFO("Warning, unexpected path length")
		return -1;
	}

	size_t maxEdgesMatched = 0;
	for (size_t i = 0; i < genomePath.size() - myPath.size() + 1; ++i) {

		size_t j = i;
		size_t edgesMatched = 0;
		size_t lengthMatched = 0;

		for(auto iter = myPath.begin(); iter != myPath.end(); ++iter) {
			if (*iter == genomePath[j]) {
				++edgesMatched;
				lengthMatched += g.length(*iter);
			}
			++j;
		}

		if (edgesMatched > maxEdgesMatched) {
			maxEdgesMatched = edgesMatched;
			maxLengthMached = lengthMatched;
			startPos = i;
		}
	}

	return maxEdgesMatched;
}


//Count all paths in genome paths
template<size_t k>
size_t PathsInGenome(Graph& g, const EdgeIndex<k + 1, Graph>& index, const Sequence& genome, std::vector<BidirectionalPath>& paths) {
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k> (genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k> (!genome, g, index);

	size_t pathCount = 0;
	for(auto iter = paths.begin(); iter != paths.end(); ++iter) {

		int s = FindInGenomePath(*iter, path1);
		if (s != -1) {
			++pathCount;
			INFO("Path of length " << PathLength(g, *iter)  << " with " << iter->size() << " edges is found in genome path starting from edge " << s)
		}

		else {
			s = FindInGenomePath(*iter, path2);
			if (s != -1) {
				++pathCount;
				INFO("Path of length " << PathLength(g, *iter) << " with " << iter->size() << " edges is found in !genome path starting from edge " << s)
			}

			else {
				int pos1 = 0, pos2 = 0;
				size_t len1 = 0, len2 = 0;
				size_t edges1 = FindInGenomeInexact(g, *iter, path1, pos1, len1);
				size_t edges2 = FindInGenomeInexact(g, *iter, path2, pos2, len2);

				if (edges1 > edges2) {
					INFO("Path partly found, percent of edges matched " << (double) edges1 / iter->size() <<
							", length percentage " << (double) len1 / PathLength(g, *iter));
				}
				else {
					INFO("Path partly found, percent of edges matched " << (double) edges2 / iter->size() <<
												", length percentage " << (double) len2 / PathLength(g, *iter));
				}

			}
		}
	}
	return pathCount;
}


} // namespace long_contigs


#endif /* QUALITY_HPP_ */
