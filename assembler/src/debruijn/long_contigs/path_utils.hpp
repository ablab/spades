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

//Increase path lengths
void IncreaseLengths(Graph& g, PathLengths& lengths, EdgeId edge, bool forward) {
	size_t len = g.length(edge);
	for(auto iter = lengths.begin(); iter != lengths.end(); ++iter) {
		*iter += len;
	}

	if (forward) {
		lengths.push_back(len);
	} else {
		lengths.push_front(0);
	}
}

//Recounting lengths form all edges to path's end
void RecountLengthsForward(Graph& g, BidirectionalPath& path, PathLengths& lengths) {
	lengths.clear();
	double currentLength = 0;

	for(auto iter = path.rbegin(); iter != path.rend(); ++iter) {
		currentLength += g.length(*iter);
		lengths.push_front(currentLength);
	}
}

//Recounting lengths from path's start to all edges
void RecountLengthsBackward(Graph& g, BidirectionalPath& path, PathLengths& lengths) {
	lengths.clear();
	double currentLength = 0;

	for(auto iter = path.begin(); iter != path.end(); ++iter) {
		lengths.push_back(currentLength);
		currentLength += g.length(*iter);
	}
}

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

template<class T>
bool ContainsPath(const T& path, const EdgeId sample) {
	for (size_t i = 0; i < path.size(); ++i) {
		if (sample == path[i]) {
			return true;
		}
	}
	return false;
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

bool ComplementPaths(Graph& g, const BidirectionalPath& path1, const BidirectionalPath& path2) {
	if (path1.size() != path2.size()) {
		return false;
	}

	for (size_t i = 0; i < path1.size(); ++i) {
		if (path1[i] != g.conjugate(path2[path1.size() - i - 1])) {
			return false;
		}
	}

	return true;
}

bool ContainsComplementPath(Graph& g, const BidirectionalPath& path, const BidirectionalPath& sample) {
	if (path.size() < sample.size()) {
		return false;
	}

	for (size_t i = 0; i < path.size() - sample.size() + 1 ; ++i) {
		bool found = true;

		for (size_t j = 0; j < sample.size(); ++j) {
			if (g.conjugate(sample[sample.size() - j - 1]) != path[i + j]) {
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

size_t LengthComplement(Graph& g, const BidirectionalPath& path, const BidirectionalPath& sample) {
	if (path.size() < sample.size()) {
		return false;
	}

	size_t max = 0;
	for (size_t i = 0; i < path.size() - sample.size() + 1 ; ++i) {
		size_t found = 0;

		for (size_t j = 0; j < sample.size(); ++j) {
			if (g.conjugate(sample[sample.size() - j - 1]) != path[i + j]) {
				continue;
			}
			found += g.length(path[i + j]);
		}

		if (max < found) {
			max = found;
		}
	}
	return max;
}

template<class T>
bool ContainsAnyOf(const BidirectionalPath& path, T& pathCollection) {
	for (auto iter = pathCollection.begin(); iter != pathCollection.end(); ++iter) {
		if (ContainsPath(path, *iter)) {
			return true;
		}
	}
	return false;
}

bool ContainsPathAt(const BidirectionalPath& path, const EdgeId sample, size_t at = 0) {
	if (at >= path.size()) {
		return false;
	}

	return path[at] == sample;
}

bool ContainsPathAt(const BidirectionalPath& path, const BidirectionalPath& sample, size_t at = 0) {
	if (path.size() < at + sample.size()) {
		return false;
	}

	for (size_t j = 0; j < sample.size(); ++j) {
		if (sample[j] != path[at + j]) {
			return false;
		}
	}

	return true;
}

template<class T>
bool ContainsAnyAt(const BidirectionalPath& path, T& pathCollection, size_t at = 0) {
	for (auto iter = pathCollection.begin(); iter != pathCollection.end(); ++iter) {
		if (ContainsPathAt(path, *iter, at)) {
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

void FilterComlementEdges(Graph& g, std::set<EdgeId>& filtered, std::set<EdgeId>& rest) {
	size_t edges = 0;
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		++edges;
		if (rest.count(*iter) == 0) {
			filtered.insert(*iter);
			if (g.conjugate(*iter) != *iter) {
				rest.insert(g.conjugate(*iter));
			}
		}
	}
	INFO("Edges separated by " << filtered.size() << " and " << rest.size() << " from " << edges);
}

void FilterComlementEdges(Graph& g, std::set<EdgeId>& filtered) {
	std::set<EdgeId> rest;
	FilterComlementEdges(g, filtered, rest);
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


void FilterUntrustedSeeds(Graph& g, std::vector<BidirectionalPath>& paths,
		std::vector<BidirectionalPath>& output, PairedInfoIndices& pairedInfo) {


}


//Remove duplicate paths
void RemoveDuplicate(Graph& g, const std::vector<BidirectionalPath>& paths,
		std::vector<BidirectionalPath>& output,
		std::vector<double>* quality = 0) {

	std::vector<BidirectionalPath> temp(paths.size());
	std::copy(paths.begin(), paths.end(), temp.begin());

	SimplePathComparator pathComparator(g);
	std::sort(temp.begin(), temp.end(), pathComparator);

	output.clear();
	if (quality != 0) {
		quality->clear();
	}

	for (int i = 0; i < (int) temp.size(); ++i) {
		bool copy = true;
		for (int j = 0; j < (int) output.size(); ++j) {
			if (ComparePaths(output[j], temp[i])) {
				copy = false;
				if (quality != 0) {
					quality->at(j) += 1.0;
				}
				break;
			}
		}

		if (copy) {
			output.push_back(temp[i]);
			if (quality != 0) {
				quality->push_back(1.0);
			}
		}
	}
}

//Remove subpaths
void RemoveSubpaths(Graph& g, std::vector<BidirectionalPath>& paths,
		std::vector<BidirectionalPath>& output,
		std::vector<double>* quality = 0) {

	std::vector<BidirectionalPath> temp(paths.size());
	std::copy(paths.begin(), paths.end(), temp.begin());

	SimplePathComparator pathComparator(g);
	std::sort(temp.begin(), temp.end(), pathComparator);
	std::vector<size_t> lengths;
	CountPathLengths(g, temp, lengths);

	output.clear();
	if (quality != 0) {
		quality->clear();
	}

	for (int i = 0; i < (int) temp.size(); ++i) {
		bool copy = true;
		for (int j = 0; j < (int) output.size(); ++j) {
			if (ContainsPath(output[j], temp[i])) {
				copy = false;
				if (quality != 0) {
					quality->at(j) += ((double) lengths[i]) / ((double) lengths[j]);
				}
				break;
			}
		}

		if (copy) {
			output.push_back(temp[i]);
			if (quality != 0) {
				quality->push_back(1.0);
			}
		}
	}
}

typedef std::multiset<EdgeId> EdgeStat;

void CountSimilarity(Graph& g, EdgeStat& path1, EdgeStat& path2, int& similarEdges, int& similarLen) {
	similarEdges = 0;
	similarLen = 0;

	auto iter = path1.begin();
	while (iter != path1.end()) {
		int count = std::min(path1.count(*iter), path2.count(*iter));
		similarEdges += count;
		similarLen += count * g.length(*iter);
		EdgeId current = *iter;
		while (iter != path1.end() && current == *iter) {
			++iter;
		}
	}

}

//Remove similar paths
void RemoveSimilar(Graph& g, std::vector<BidirectionalPath>& paths,
		std::vector<double>& quality,
		std::set<int>& toRemove) {

	INFO("Removing similar");

	toRemove.clear();
	std::vector<EdgeStat> pathStat;

	std::vector<size_t> lengths;
	CountPathLengths(g, paths, lengths);

	DETAILED_INFO("Counting stats");
	for (int i = 0; i < (int) paths.size(); ++i) {
		EdgeStat stat;
		for (auto edge = paths[i].begin(); edge != paths[i].end(); ++edge) {
			stat.insert(*edge);
		}
		pathStat.push_back(stat);
	}
	for (int i = 0; i < (int) paths.size(); ++i) {
		if (toRemove.count(i) != 0) {
			continue;
		}

		for (int j = i + 1; j < (int) paths.size(); ++j) {
			int similarLen = 0;
			int similarEdges = 0;
			CountSimilarity(g, pathStat[i], pathStat[j], similarEdges, similarLen);

			if (((double) similarLen) / ((double) lengths[j]) >= lc_cfg::get().fo.similar_length &&
					((double) similarEdges) / ((double) paths[j].size()) >= lc_cfg::get().fo.similar_edges) {
				toRemove.insert(j);
			}
		}
	}

	INFO("Done");
}


bool HasConjugate(Graph& g, BidirectionalPath& path) {
	size_t count = 0;
	double len = 0.0;

	for (auto e1 = path.begin(); e1 != path.end(); ++e1) {
		for (auto e2 = path.begin(); e2 != path.end(); ++e2) {
			if (*e1 != *e2 && g.conjugate(*e1) == *e2) {
				++count;
				len += g.length(*e2);
				break;
			}
		}
	}

	if (count != 0) {
		DETAILED_INFO("Self conjugate detected: edges " << count << ", length: " << len);
		DetailedPrintPath(g, path);
	}

	return math::gr(len / PathLength(g, path), lc_cfg::get().fo.conj_len_percent);
}


void BreakApart(Graph& g, BidirectionalPath& path, std::vector<BidirectionalPath>& output) {
	int i, j;
	bool found = false;

	for (i = 0; i < (int) path.size(); ++i) {
		for (j = path.size() - 1; j > i; --j) {
			if (g.conjugate(path[i]) == path[j]) {
				found = true;
				break;
			}
		}
		if (found) {
			break;
		}
	}

	if (!found) {
		return;
	}

	while (g.conjugate(path[i]) == path[j]) {
		++i;
		--j;
	}

	for (int k = i; k <= j; ++k) {
		if (g.length(path[k]) >= K - lc_cfg::get().fo.chimeric_delta &&
				g.length(path[k]) <= K + lc_cfg::get().fo.chimeric_delta) {

			BidirectionalPath left;
			for (int l = 0; l < k; ++l) {
				left.push_back(path[l]);
			}
			BidirectionalPath right;
			for (int l = k + 1; l < (int) path.size(); ++l) {
				right.push_back(path[l]);
			}

			output.push_back(left);
			output.push_back(right);

			return;
		}
	}

	BidirectionalPath left;
	for (int l = 0; l < i; ++l) {
		left.push_back(path[l]);
	}
	BidirectionalPath right;
	for (int l = j + 1; l < (int) path.size(); ++l) {
		right.push_back(path[l]);
	}

	output.push_back(left);
	output.push_back(right);
}

void RemoveWrongConjugatePaths(Graph& g, std::vector<BidirectionalPath>& paths,
		std::vector<BidirectionalPath>& output) {

	output.clear();
	for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
		if (!HasConjugate(g, *iter)) {
			output.push_back(*iter);
		} else {
			INFO("Removed as self conjugate");
			if (lc_cfg::get().fo.break_sc) {
				INFO("Added half");
				BreakApart(g, *iter, output);
			}
		}
	}
}


std::pair<int, int> FindSameSearchRange(std::vector<BidirectionalPath>& paths, std::vector<size_t>& lengths, int pathNum) {
	size_t length = lengths[pathNum];

	int i = pathNum - 1;
	while (i >= 0 && lengths[i] == length) {
		--i;
	}

	int j = pathNum + 1;
	while (j < (int) paths.size() && lengths[j] == length) {
		++j;
	}

	return std::make_pair(i + 1, j - 1);
}


std::pair<int, int> FindSearchRange(std::vector<BidirectionalPath>& paths, std::vector<size_t>& lengths, int pathNum) {
	static double coeff = lc_cfg::get().fo.length_percent;
	double length = (double) lengths[pathNum];

	int i = pathNum - 1;
	while (i >=0 && lengths[i] <= length * coeff) {
		--i;
	}

	int j = pathNum + 1;
	while (j < (int) paths.size() && lengths[j] >= length / coeff) {
		++j;
	}

	return std::make_pair(i + 1, j - 1);
}


std::pair<int, double> FindComlementPath(Graph& g, std::vector<BidirectionalPath>& paths, std::vector<size_t>& lengths, int pathNum) {
	static double conjugatePercent = lc_cfg::get().fo.conjugate_percent;
	BidirectionalPath& path = paths[pathNum];

	auto range = FindSameSearchRange(paths, lengths, pathNum);
	for (int i = range.first; i <= range.second; ++i) {
		if (ComplementPaths(g, path, paths[i])) {
			INFO("Total complemented");
			return std::make_pair(i, 1.0);
		}
	}

	double maxConj = 0;
	int maxI = 0;

	range = FindSearchRange(paths, lengths, pathNum);
	for (int i = range.first; i <= range.second; ++i) {
		BidirectionalPath& p = (path.size() < paths[i].size()) ? paths[i] : path;
		BidirectionalPath& sample = (path.size() < paths[i].size()) ? path : paths[i];

		if (ContainsComplementPath(g, p, sample)) {
			INFO("Complement subpath " << ((double) PathLength(g, sample)) / ((double) PathLength(g, p)));
			PrintPath(g, p);
			PrintPath(g, sample);
			return std::make_pair(i, ((double) PathLength(g, sample)) / ((double) PathLength(g, p)));
		}

		size_t conjLength = LengthComplement(g, p, sample);
		double c = ((double) conjLength) / ((double) PathLength(g, p));

		if (c > maxConj) {
			maxConj = c;
			maxI = i;
		}
	}

	if (maxConj >= conjugatePercent) {
		INFO("Partly complement with " << maxI << ", percentage " << maxConj);
		PrintPath(g, path);
		PrintPath(g, paths[maxI]);
		return std::make_pair(maxI, maxConj);
	}
	else {
		INFO("NO COMPLEMENT!");
		PrintPath(g, path);
		return std::make_pair(-1, 0);
	}
}

//Filter symetric complement contigs
void FilterComplement(Graph& g, std::vector<BidirectionalPath>& paths, std::vector<int>* pairs, std::vector<double>* quality) {

	std::sort(paths.begin(), paths.end(), SimplePathComparator(g));
	pairs->clear();
	quality->clear();
	pairs->resize(paths.size(), -1);
	quality->resize(paths.size(), 0);

	std::vector<size_t> lengths;
	CountPathLengths(g, paths, lengths);

	std::set<int> found;

	int i = 0;
	int revert = -1;
	while (i < (int) paths.size()) {
		if (found.count(i) == 0) {
			auto comp = FindComlementPath(g, paths, lengths, i);

			if (found.count(comp.first) > 0 && comp.first != i && comp.first != -1) {
				INFO("Wrong complement pairing");
				PrintPath(g, paths[i]);
				PrintPath(g, paths[pairs->at(comp.first)]);
				PrintPath(g, paths[comp.first]);
				INFO("Substituting");

				if (comp.second <= quality->at(comp.first)) {
					++i;
					continue;
				}
				INFO("Will substitute. New quality " << comp.second << " greater than " << quality->at(comp.first));
				found.erase(pairs->at(comp.first));
				pairs->at(pairs->at(comp.first)) = -1;
				quality->at(pairs->at(comp.first)) = 0;
				revert = pairs->at(comp.first);
			}

			if (comp.first == -1) {
				INFO("Really not found");
				++i;
				continue;
			}

			found.insert(i);
			found.insert(comp.first);
			pairs->at(i) = comp.first;
			pairs->at(comp.first) = i;
			quality->at(i) = comp.second;
			quality->at(comp.first) = comp.second;
		}
		++i;

		if (revert != -1 && i > revert) {
			INFO("Reverting from " << i << " to " << revert);
			i = revert;
			revert = -1;
		}
	}

	INFO("Results of complement filtering")
	found.clear();
	for (int i = 0; i < (int) paths.size(); ++i) {
		if (found.count(i) == 0) {
			if (quality->at(i) == 1) {
				INFO("Total complement");
			}
			else {
				INFO("Complement subpath " << quality->at(i));
				PrintPath(g, paths[i]);
				if (pairs->at(i) != -1) {
					PrintPath(g, paths[pairs->at(i)]);
				} else {
					INFO("No complment path");
				}
			}
			found.insert(i);
			found.insert(pairs->at(i));
		}
	}
}

//Remove overlaps, remove sub paths first
void RemoveOverlaps(Graph& g, std::vector<BidirectionalPath>& paths, std::vector<int>& pairs, std::vector<double>& quality) {
	INFO("Removing overlaps");
    for (int k = 0; k < (int) paths.size(); ++k) {
    		BidirectionalPath& path = paths[k];
            EdgeId lastEdge = path.back();

            int overlap = -1;
            int overlaped = - 1;
            for (int l = 0; l < (int) paths.size(); ++l) {
				if (k != l) {
					BidirectionalPath& toCompare = paths[l];

					for (int i = 0; i < (int) toCompare.size(); ++i) {
						if (lastEdge == toCompare[i]) {
							int diff = path.size() - i - 1;
							bool found = true;

							for (int j = i - 1; j >= 0; --j) {
								if (toCompare[j] != path[j + diff]) {
									found = false;
									break;
								}
							}

							if (found && overlap < i) {
								overlap = i;
								overlaped = l;
							}
					   }
					}
				}
            }

            if (overlap != -1) {
            	size_t overlapLength = 0;
				for (int i = 0; i <= overlap; ++i) {
					overlapLength += g.length(path.back());
				}

            	INFO("Found overlap by " << overlap + 1 << " edge(s) with total length " << overlapLength);
            	PrintPath(g, path);
            	PrintPath(g, paths[overlaped]);

				overlap = std::min(overlap, (int) path.size() - 1);

				for (int i = 0; i <= overlap; ++i) {
					path.pop_back();
				}

				if (quality[k] == 1.0) {
					INFO("Same one removed from reverse-complement path");
					BidirectionalPath& comp = paths[pairs[k]];
					for (int i = 0; i <= overlap; ++i) {
						comp.pop_front();
					}
				}
            }
    }
	INFO("Done");
}

void MakeBlackSet(Graph& g, Path<Graph::EdgeId>& path1, Path<Graph::EdgeId>& path2, std::set<EdgeId> blackSet) {
	for (auto edge1 = g.SmartEdgeBegin(); !edge1.IsEnd(); ++edge1) {
		if (!ContainsPath(path1, *edge1) && !ContainsPath(path2, *edge1)) {
			blackSet.insert(*edge1);
		}
	}
}


} // namespace long_contigs

#endif /* PATH_UTILS_HPP_ */
