//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * weight_counter.hpp
 *
 *  Created on: Feb 19, 2012
 *      Author: andrey
 */

#ifndef WEIGHT_COUNTER_HPP_
#define WEIGHT_COUNTER_HPP_

#include "bidirectional_path.hpp"
#include "paired_library.hpp"
#include <algorithm>
#include <boost/math/special_functions/fpclassify.hpp>

namespace path_extend {

struct EdgeWithPairedInfo {
	size_t e_;
	double pi_;

	EdgeWithPairedInfo(size_t e_, double pi) :
			e_(e_), pi_(pi) {

	}
};

struct EdgeWithDistance {
	EdgeId e_;
	int d_;

	EdgeWithDistance(EdgeId e, size_t d) :
			e_(e), d_(d) {
	}
};

class ExtentionAnalyzer {

protected:
	Graph& g_;
	PairedInfoLibrary& lib_;

public:
	ExtentionAnalyzer(Graph& g, PairedInfoLibrary& l) :
			g_(g), lib_(l) {
	}

	PairedInfoLibrary& getLib() {
		return lib_;
	}

	void FindCoveredEdges(const BidirectionalPath& path, EdgeId candidate,
			std::vector<EdgeWithPairedInfo>& edges) {
		edges.clear();

		for (int i = (int) path.Size() - 1; i >= 0; --i) {
			double w = lib_.IdealPairedInfo(path[i], candidate,
					path.LengthAt(i));
			if (math::gr(w, 0.)) {
			}
				edges.push_back(EdgeWithPairedInfo(i, w));
		}
	}

	void FindForwardEdges(const BidirectionalPath& path, EdgeId candidate,
			std::vector<EdgeWithDistance>& edges) {
		edges.clear();
		edges.push_back(EdgeWithDistance(candidate, 0));

		size_t i = 0;
		while (i < edges.size()) {
			size_t currentDistance = edges[i].d_ + g_.length(edges[i].e_);
			auto nextEdges = g_.OutgoingEdges(g_.EdgeEnd(edges[i].e_));

			if (edges[i].d_ < (int) lib_.insert_size_) {
				for (auto edge = nextEdges.begin(); edge != nextEdges.end();
						++edge) {
					edges.push_back(EdgeWithDistance(*edge, currentDistance));
				}
			}
			++i;
		}
	}
};

class WeightCounter {

protected:
	Graph& g_;
	PairedInfoLibraries& libs_;
	std::vector<ExtentionAnalyzer *> analyzers_;
	double avrageLibWeight_;

	double threshold_;
	bool normalizeWeight_;
	bool normalizeWightByCoverage_;

	bool tryDeepSearch_;

	std::set<int> excludedEdges_;

public:

	WeightCounter(Graph& g, PairedInfoLibraries& libs, double threshold = 0.0) :
			g_(g), libs_(libs), threshold_(threshold), normalizeWeight_(true), normalizeWightByCoverage_(
					true), tryDeepSearch_(false), excludedEdges_() {
		avrageLibWeight_ = 0.0;
		analyzers_.reserve(libs_.size());
		for (auto iter = libs_.begin(); iter != libs_.end(); ++iter) {
			analyzers_.push_back(new ExtentionAnalyzer(g_, **iter));
			avrageLibWeight_ += (*iter)->coverage_coeff_;
		}
		avrageLibWeight_ /= max(libs_.size(), (size_t) 1);
	}

	virtual ~WeightCounter() {
		for (auto iter = analyzers_.begin(); iter != analyzers_.end(); ++iter) {
			delete *iter;
		}
		analyzers_.clear();
	}

	virtual bool PairInfoExist(EdgeId first, EdgeId second, int distance) = 0;

	virtual double CountWeight(BidirectionalPath& path, EdgeId e,
			int gapLength = 0) = 0;

	virtual void GetDistances(EdgeId e1, EdgeId e2, std::vector<int>& dist,
			std::vector<double>& w) = 0;

	virtual double CountIdealInfo(EdgeId e1, EdgeId e2, size_t dist) = 0;

	virtual double CountIdealInfo(const BidirectionalPath& p, EdgeId e,
			size_t gap) = 0;

	virtual bool IsExtensionPossible(BidirectionalPath& path, EdgeId e) {
		return IsExtensionPossible(CountWeight(path, e)) ? true : false;
	}

	virtual bool IsExtensionPossible(double weight) {
		return math::ge(weight, threshold_) ? true : false;
	}

	std::set<int> & GetExcludedEdges() {
		return excludedEdges_;
	}

	double getThreshold() const {
		return threshold_;
	}

	bool isNormalizeWeight() const {
		return normalizeWeight_;
	}

	bool isNormalizeWightByCoverage() const {
		return normalizeWightByCoverage_;
	}

	void setNormalizeWeight(bool normalizeWeight) {
		this->normalizeWeight_ = normalizeWeight;
	}

	void setNormalizeWightByCoverage(bool normalizeWightByCoverage) {
		this->normalizeWightByCoverage_ = normalizeWightByCoverage;
	}

	void setThreshold(double threshold) {
		this->threshold_ = threshold;
	}

	void setDeepSearch(bool value) {
		tryDeepSearch_ = value;
	}

	bool isDeepSearch() {
		return tryDeepSearch_;
	}

	PairedInfoLibraries& getLibs() {
		return libs_;
	}

};

class ReadCountWeightCounter: public WeightCounter {

protected:

	double CountSingleLib(int libIndex, BidirectionalPath& path, EdgeId e,
			int additionalGapLength = 0.0) {

		PairedInfoLibrary& pairedInfoLibrary = *libs_[libIndex];
		double weight = 0.0;

		std::vector<EdgeWithPairedInfo> coveredEdges;
		analyzers_[libIndex]->FindCoveredEdges(path, e, coveredEdges);

		for (auto iter = coveredEdges.begin(); iter != coveredEdges.end();
				++iter) {
			if (excludedEdges_.count(iter->e_) > 0) {
				DEBUG("excluded " << iter->e_)
				continue;
			}
			double w = libs_[libIndex]->CountPairedInfo(path[iter->e_], e,
					path.LengthAt(iter->e_) + additionalGapLength);

			if (normalizeWeight_) {
				w /= iter->pi_;
			}
			weight += w;

		}

		return normalizeWightByCoverage_ ?
				pairedInfoLibrary.normalizeByCoverage(weight)
						* avrageLibWeight_ :
				weight;
	}

public:

	ReadCountWeightCounter(Graph& g_, PairedInfoLibraries& libs_,
			double threshold_ = 0.0) :
			WeightCounter(g_, libs_, threshold_) {
	}

	virtual void GetDistances(EdgeId e1, EdgeId e2, std::vector<int>& dist,
			std::vector<double>& w) {
		for (size_t i = 0; i < libs_.size(); ++i) {
			libs_[i]->CountDistances(e1, e2, dist, w);
		}
	}

	virtual double CountIdealInfo(EdgeId e1, EdgeId e2, size_t dist) {
		double res = 0.0;
		for (size_t i = 0; i < libs_.size(); ++i) {
			res += libs_[i]->IdealPairedInfo(e1, e2, dist);
		}
		return res;
	}

	virtual double CountIdealInfo(const BidirectionalPath& p, EdgeId e,
			size_t gap) {
		double w = 0.0;
		for (int i = (int) p.Size() - 1; i >= 0; --i) {
			w += CountIdealInfo(p[i], e, gap + p.LengthAt(i));
		}
		return w;
	}

	virtual double CountWeight(BidirectionalPath& path, EdgeId e,
			int gapLength = 0) {
		double weight = 0.0;
		std::vector<EdgeWithDistance> edges;

		for (size_t i = 0; i < libs_.size(); ++i) {
			if (tryDeepSearch_) {
				analyzers_[i]->FindForwardEdges(path, e, edges);

				for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
					weight += CountSingleLib(i, path, edge->e_,
							edge->d_ + gapLength);
				}
			} else {
				weight += CountSingleLib(i, path, e, gapLength);
			}
		}

		return weight;
	}

	virtual bool PairInfoExist(EdgeId first, EdgeId second, int distance) {
		for (size_t libIndex = 0; libIndex < libs_.size(); ++libIndex) {
			double w = libs_[libIndex]->CountPairedInfo(first, second,
					distance);
			double w_ideal = libs_[libIndex]->IdealPairedInfo(first, second,
					distance);
			if (w_ideal == 0) {
				continue;
			}
			if (normalizeWeight_) {
				w /= w_ideal;
			}
			if (w > 0) {
				return true;
			}
		}
		return false;
	}

};

class PathCoverWeightCounter: public WeightCounter {

protected:

	double singleThreshold;

	double CountSingleLib(int libIndex, BidirectionalPath& path, EdgeId e,
			int additionalGapLength = 0.0) {
		PairedInfoLibrary& pairedInfoLibrary = *libs_[libIndex];
		double weight = 0.0;
		double idealWeight = 0.0;

		std::vector<EdgeWithPairedInfo> coveredEdges;

		analyzers_[libIndex]->FindCoveredEdges(path, e, coveredEdges);
		std::set<int> excludedEdges(excludedEdges_);
		for (auto iter = coveredEdges.begin(); iter != coveredEdges.end();
				++iter) {
			if (excludedEdges.count(iter->e_) > 0) {
				excludedEdges.erase(excludedEdges.find(iter->e_));
				continue;
			}

			double threshold =
					pairedInfoLibrary.single_threshold_ >= 0.0 ?
							pairedInfoLibrary.single_threshold_ :
							singleThreshold;

			double singleWeight = libs_[libIndex]->CountPairedInfo(
					path[iter->e_], e,
					path.LengthAt(iter->e_) + additionalGapLength);

			if (normalizeWeight_) {
				singleWeight /= iter->pi_;
			}
			if (normalizeWightByCoverage_) {
				singleWeight = pairedInfoLibrary.normalizeByCoverage(
						singleWeight) * avrageLibWeight_;
			}
			if (math::ge(singleWeight, threshold)) {
				weight += iter->pi_;
			}
			idealWeight += iter->pi_;
		}

		return math::gr(idealWeight, 0.0) ? weight / idealWeight : 0.0;
	}


public:

	PathCoverWeightCounter(Graph& g_, PairedInfoLibraries& libs_,
			double threshold_ = 0.0, double singleThreshold_ = 0.0) :
			WeightCounter(g_, libs_, threshold_), singleThreshold(
					singleThreshold_) {

	}

	virtual void GetDistances(EdgeId e1, EdgeId e2, std::vector<int>& dist,
			std::vector<double>& w) {
		for (size_t i = 0; i < libs_.size(); ++i) {
			libs_[i]->CountDistances(e1, e2, dist, w);
		}
	}

	virtual double CountIdealInfo(EdgeId e1, EdgeId e2, size_t dist) {
		double res = 0.0;
		for (size_t i = 0; i < libs_.size(); ++i) {
			res += libs_[i]->IdealPairedInfo(e1, e2, dist);
		}
		return res;
	}

	virtual double CountIdealInfo(const BidirectionalPath& p, EdgeId e,
			size_t gap) {
		double w = 0.0;
		for (int i = (int) p.Size() - 1; i >= 0; --i) {
			w += g_.length(p[i]) ?
					CountIdealInfo(p[i], e, gap + p.LengthAt(i)) > 0 : 0;
		}
		return w;
	}

	virtual double CountWeight(BidirectionalPath& path, EdgeId e,
			int gapLength = 0) {
		double weight = 0.0;
		for (size_t i = 0; i < libs_.size(); ++i) {
			weight += CountSingleLib(i, path, e, gapLength);
		}

		return weight / max(libs_.size(), (size_t) 1);
	}

	virtual bool PairInfoExist(EdgeId first, EdgeId second, int distance) {
		for (size_t libIndex = 0; libIndex < libs_.size(); ++libIndex) {
			double w = libs_[libIndex]->CountPairedInfo(first, second,
					distance);
			double w_ideal = libs_[libIndex]->IdealPairedInfo(first, second,
					distance);
			if (w_ideal == 0.0) {
				continue;
			}
			if (normalizeWeight_) {
				w /= w_ideal;
			}
			double threshold =
					libs_[libIndex]->single_threshold_ >= 0.0 ?
							libs_[libIndex]->single_threshold_ :
							singleThreshold;
			DEBUG("pair info exsist " <<g_.int_id(first) << "  " << g_.int_id(second) << " " << w << " " << distance << " " << threshold << " " << w_ideal);
			if (w > threshold) {
				return true;
			}
		}
		return false;
	}

};

class LongReadsWeightCounter: public WeightCounter {

    size_t RL_;
    GraphCoverageMap& coverageMap_;


public:

    LongReadsWeightCounter(Graph& g_,  GraphCoverageMap& covMap, size_t RL,
            double threshold = 0.0): WeightCounter(g_, *(new std::vector<PairedInfoLibrary *> ()), threshold), RL_(RL), coverageMap_(covMap)
            {

    }

    virtual void GetDistances(EdgeId e1, EdgeId e2, std::vector<int>& dist,
            std::vector<double>& w) {
    }

    virtual double CountIdealInfo(EdgeId e1, EdgeId e2, size_t dist) {
        int res = std::min(g_.length(e1) + dist + g_.length(e2) - RL_, g_.length(e1) - 1) - std::max(g_.length(e1) + dist + 1 - RL_, (size_t)0) + 1;
        return res > 0 ? res : 0;
    }

    virtual double CountIdealInfo(const BidirectionalPath& p, EdgeId e,
            size_t gap) {
        double w = 0.0;
        for (int i = (int) p.Size() - 1; i >= 0; --i) {
            w += g_.length(p[i]) ? CountIdealInfo(p[i], e, gap + p.LengthAt(i)) > 0 : 0;
        }
        return w;
    }

    virtual double CountWeight(BidirectionalPath& path, EdgeId e,
            int gapLength = 0)
    {
        std::vector<EdgeWithPairedInfo> coveredEdges;
        std::set<int> includedEdges;
        std::set<int> allEdges;
        for (int index = (int) path.Size() - 1; index >=0; --index) {
        	allEdges.insert(index);
        	if (excludedEdges_.count(index) > 0) {
            	continue;
            }
        	EdgeId cur_edge = path.At(index);
        	double idealWeight = IdealWeightBetweenEdges(cur_edge, e, path.LengthAt(index));
        	if (idealWeight > 0){
        		double singleWeight = Weight(cur_edge, e, path.LengthAt(index));
        		singleWeight /= idealWeight;
        		if (math::ge(singleWeight, threshold_)){
        			includedEdges.insert(index);
        		}
        	}
        }
        double weight = CountIdealWeight(path, e, includedEdges);
        double commonIdealWeight = CountIdealWeight(path, e, allEdges);

        return math::gr(commonIdealWeight, 0.0) ? weight / commonIdealWeight : 0.0;
    }

    double CountIdealWeight(BidirectionalPath& p, EdgeId e1,
			set<int>& includedEdges) {
    	size_t common_length = 0;
    	double ideal_length = RL_ - g_.k() - 1;
    	double weight =  min(ideal_length, (double)p.Length() - 1.);
		for (int index = (int) p.Size() - 1; index >= 0; --index) {
			EdgeId cur_edge = p.At(index);
			if (common_length > RL_ - g_.k() - 1){
				break;
			}
			if (includedEdges.count(index) > 0) {
				common_length =  p.LengthAt(index);
				continue;
			}
			common_length -= min(g_.length(cur_edge), RL_ - g_.k() - 1 - common_length);
			common_length =  p.LengthAt(index);
		}
		return weight;
	}

    double IdealWeightBetweenEdges(EdgeId e1, EdgeId e2, size_t dist){
    	double w = RL_ - g_.k() - dist - 2 + g_.length(e1);
    	return std::max(w, 0.);
    }

    double Weight(EdgeId e1, EdgeId e2, size_t gap) {
           auto supp_paths = coverageMap_.GetCoveringPaths(e1);
           double weight = 0;
           for (auto it = supp_paths.begin(); it != supp_paths.end(); ++it) {
               BidirectionalPath* cur_path = *it;
               auto poss = cur_path->FindAll(e1);
               for (size_t i = 0; i < poss.size(); ++i) {
                   size_t pos = poss[i];
                   size_t cur_gap = 0;
                   while (pos < cur_path->Size() && cur_gap <= gap){
                       EdgeId cur_edge = cur_path->At(pos);
                       if (cur_edge == e2 && cur_gap == gap){
                           weight += cur_path->getWeight();
                           break;
                       }
                       cur_gap += g_.length(cur_edge);
                       pos++;
                   }
               }
           }
           return weight;
       }

    virtual bool PairInfoExist(EdgeId first, EdgeId second, int distance) {
        return false;
    }

};

}

#endif /* WEIGHT_COUNTER_HPP_ */
