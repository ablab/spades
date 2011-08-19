/*
 * path_utils.hpp
 *
 *  Created on: Aug 17, 2011
 *      Author: andrey
 */

#ifndef PATH_UTILS_HPP_
#define PATH_UTILS_HPP_

#include "lc_common.hpp"

namespace long_contigs {

bool ComparePaths(const BidirectionalPath& path1, const BidirectionalPath& path2) {
	if (path1.size() != path2.size()) {
		return false;
	}

	for (size_t i = 0; i < path1.size(); ++i) {
		if (path1[i] != path2[i]) {
			return false;
		}
	}
	return true;
}

bool ContainsPath(const BidirectionalPath& path, const BidirectionalPath& sample) {
	if (path.size() < sample.size()) {
		return false;
	}

	for (size_t i = 0; i < path.size() - sample.size() + 1 ; ++i) {
		bool found = true;

		for (size_t j = 0; j < sample.size(); ++j) {
			if (sample[j] != path[i + j]) {
				found = false;
				break;
			}
		}

		if (found) {
			return true;
		}
	}
	return false;
}

//Find coverage of worst covered edge
double PathMinReadCoverage(Graph& g, BidirectionalPath& path) {
	if (path.empty()) {
		return 0;
	}

	double minCov = g.coverage(path[0]);

	for (auto edge = path.begin(); edge != path.end(); ++edge) {
		double cov = g.coverage(*edge);
		if (minCov > cov) {
			minCov = cov;
		}
	}

	return minCov;
}

//Filter paths by containing paths
void FilterPaths(Graph& g, std::vector<BidirectionalPath>& paths, std::vector<BidirectionalPath>& samples) {
	for(auto path = paths.begin(); path != paths.end(); ) {
		bool toErase = true;
		for (auto sample = samples.begin(); sample != samples.end(); ++sample) {
			if (ContainsPath(*path, *sample)) {
				toErase = false;
				break;
			}
		}
		if (toErase) {
			paths.erase(path);
		} else {
			++path;
		}
	}
}

//Filter paths only with edges with given length
void FilterEdge(Graph& g, std::vector<BidirectionalPath>& paths, size_t edgeLen) {
	std::vector<BidirectionalPath> samples;
	for (auto edge = g.SmartEdgeBegin(); !edge.IsEnd(); ++edge) {
		if (g.length(*edge) == edgeLen) {
			BidirectionalPath sample;
			sample.push_back(*edge);
			samples.push_back(sample);
		}
	}
	FilterPaths(g, paths, samples);
}

//Remove paths with low covered edges
void FilterLowCovered(Graph& g, std::vector<BidirectionalPath>& paths, double threshold) {
	for (auto path = paths.begin(); path != paths.end(); ) {
		if (PathMinReadCoverage(g, *path) < threshold) {
			paths.erase(path);
		} else {
			++path;
		}
	}
}


//Remove duplicate paths
void RemoveDuplicate(const std::vector<BidirectionalPath>& paths, std::vector<BidirectionalPath>& output) {
	for (auto path = paths.begin(); path != paths.end(); ++path) {
		bool copy = true;
		for (auto iter = output.begin(); iter != output.end(); ++iter) {
			if (ComparePaths(*path, *iter)) {
					copy = false;
					break;
			}
		}

		if (copy) {
			output.push_back(*path);
		}
	}
}

class SimplePathComparator {
private:
	Graph& g_;

public:
	SimplePathComparator(Graph& g): g_(g) {}

	bool operator() (const BidirectionalPath& path1, const BidirectionalPath& path2) {
		return PathLength(g_, path1) > PathLength(g_, path2);
	}
};

//Remove subpaths
void RemoveSubpaths(Graph& g, std::vector<BidirectionalPath>& paths, std::vector<BidirectionalPath>& output) {
	std::vector<BidirectionalPath> temp(paths.size());
	std::copy(paths.begin(), paths.end(), temp.begin());

	SimplePathComparator pathComparator(g);
	std::sort(temp.begin(), temp.end(), pathComparator);

	for (auto path = temp.begin(); path != temp.end(); ++path) {
		bool copy = true;
		for (auto iter = output.begin(); iter != output.end(); ++iter) {
			if (ContainsPath(*iter, *path)) {
				copy = false;
				break;
			}
		}

		if (copy) {
			output.push_back(*path);
		}
	}
}

//Remove overlaps, remove sub paths first
//TODO
void RemoveOverlaps(std::vector<BidirectionalPath>& paths) {
	INFO("Removing overlaps");
	for (auto path = paths.begin(); path != paths.end(); ++path) {
		EdgeId lastEdge = path->back();

		for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
			if (iter != path) {
				BidirectionalPath& toCompare = *iter;
				int overlap = -1;

				for (size_t i = 0; i < toCompare.size(); ++i) {
					if (lastEdge == toCompare[i]) {
						int diff = path->size() - i;
						bool found = true;

						for (int j = i - 1; j >= 0; --j) {
							if (toCompare[j] != path->at(j + diff)) {
								found = false;
								break;
							}
						}

						if (found) {
							overlap = i;
							INFO("Found overlap by " << i);
						}
					}
				}

				for (int i = 0; i <= overlap; ++i) {
					toCompare.pop_front();
				}
			}
		}
	}
	INFO("Done");
}

} // namespace long_contigs

#endif /* PATH_UTILS_HPP_ */
