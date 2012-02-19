/*
 * extension.hpp
 *
 *  Created on: Mar 5, 2012
 *      Author: andrey
 */

#ifndef EXTENSION_HPP_
#define EXTENSION_HPP_

#include "weight_counter.hpp"
#include <iostream>
#include <fstream>

namespace path_extend {

typedef std::multimap<double, EdgeWithDistance> AlternativeConteiner;


class PathAnalyzer {

protected:
    Graph& g_;

public:
    PathAnalyzer(Graph& g): g_(g) {
    }

    int ExcludeTrivial(const BidirectionalPath& path, std::set<int>& edges, int from = -1) {
        edges.clear();

        int edgeIndex = (from == -1) ? path.Size() - 1 : from;
        if ((int) path.Size() <= from) {
            return edgeIndex;
        }

        VertexId currentVertex = g_.EdgeEnd(path[edgeIndex]);
        while (edgeIndex >= 0 && g_.CheckUniqueIncomingEdge(currentVertex)) {
            EdgeId e = g_.GetUniqueIncomingEdge(currentVertex);
            currentVertex = g_.EdgeStart(e);

            edges.insert(edgeIndex);
            --edgeIndex;
        }

        return edgeIndex;
    }

    int ExcludeTrivialWithBulges(const BidirectionalPath& path, std::set<int>& edges) {
        edges.clear();

        if (path.Empty()) {
            return 0;
        }

        int lastEdge = path.Size() - 1;
        do {
            lastEdge = ExcludeTrivial(path, edges, lastEdge);

            if (lastEdge >= 0) {
                VertexId v = g_.EdgeEnd(path[lastEdge]);
                VertexId u = g_.EdgeStart(path[lastEdge]);
                auto bulgeCandidates = g_.IncomingEdges(v);
                bool bulge = true;

                for (auto iter = bulgeCandidates.begin(); iter != bulgeCandidates.end(); ++iter) {
                    if (g_.EdgeStart(*iter) != u) {
                        bulge = false;
                        break;
                    }
                }

                if (!bulge) {
                    break;
                }

                --lastEdge;
            }
        } while (lastEdge >= 0);

        return lastEdge;
    }
};


class ExtensionChooserListener {

public:

    virtual void ExtensionChosen(double weight) = 0;

    virtual void ExtensionChosen(AlternativeConteiner& alts) = 0;

    virtual ~ExtensionChooserListener() {

    }
};


class ExtensionChooser {

public:
    typedef std::vector<EdgeWithDistance> EdgeContainer;

protected:
    Graph& g_;

    WeightCounter * wc_;

    PathAnalyzer analyzer_;

    double priorityCoefficient_;

    bool excludeTrivial_;
    bool excludeTrivialWithBulges_;

    std::vector<ExtensionChooserListener *> listeners_;

public:
    ExtensionChooser(Graph& g_, WeightCounter * wc, double priority): g_(g_), wc_(wc), analyzer_(g_), priorityCoefficient_(priority),
        excludeTrivial_(true), excludeTrivialWithBulges_(true), listeners_() {
    }

    virtual ~ExtensionChooser() {

    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) = 0;

    bool isExcludeTrivial() const
    {
        return excludeTrivial_;
    }

    double CountWeight(BidirectionalPath& path, EdgeId e) {
        return wc_->CountWeight(path, e);
    }

    bool isExcludeTrivialWithBulges() const
    {
        return excludeTrivialWithBulges_;
    }

    void setExcludeTrivial(bool excludeTrivial)
    {
        this->excludeTrivial_ = excludeTrivial;
    }

    void setExcludeTrivialWithBulges(bool excludeTrivialWithBulges)
    {
        this->excludeTrivialWithBulges_ = excludeTrivialWithBulges;
    }

    PairedInfoLibraries& getLibs() {
        return wc_->getLibs();
    }

    void Subscribe(ExtensionChooserListener * listener) {
        listeners_.push_back(listener);
    }

    void NotifyAll(double weight) {
        for (auto iter = listeners_.begin(); iter != listeners_.end(); ++iter) {
            (*iter)->ExtensionChosen(weight);
        }
    }

    void NotifyAll(AlternativeConteiner& alts) {
        for (auto iter = listeners_.begin(); iter != listeners_.end(); ++iter) {
            (*iter)->ExtensionChosen(alts);
        }
    }

};


class TrivialExtensionChooser: public ExtensionChooser {

public:
    TrivialExtensionChooser(Graph& g_): ExtensionChooser(g_, 0 ,0.0)  {
    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.size() == 1) {
            return edges;
        }
        return EdgeContainer();
    }
};


class TrivialExtensionChooserWithPI: public ExtensionChooser {

public:
    TrivialExtensionChooserWithPI(Graph& g_, WeightCounter * wc): ExtensionChooser(g_, wc, 0.0) {
    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        wc_->GetExcludedEdges().clear();
        if (edges.size() == 1) {
                double weight = wc_->CountWeight(path, edges.back().e_);
                NotifyAll(weight);

                if (wc_->IsExtensionPossible(weight)) {
                    return edges;
                }
        }
        return EdgeContainer();
    }
};


class SimpleExtensionChooser: public ExtensionChooser {

protected:

    void RemoveTrivial(BidirectionalPath& path) {
        wc_->GetExcludedEdges().clear();

        if (excludeTrivialWithBulges_) {
            analyzer_.ExcludeTrivialWithBulges(path, wc_->GetExcludedEdges());
        }
        else if (excludeTrivial_) {
            analyzer_.ExcludeTrivial(path, wc_->GetExcludedEdges());
        }
    }

public:
    SimpleExtensionChooser(Graph& g, WeightCounter * wc, double priority): ExtensionChooser(g, wc, priority) {

    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.empty()) {
            return edges;
        }

        RemoveTrivial(path);
        AlternativeConteiner weights;

        for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
            double weight = wc_->CountWeight(path, iter->e_);
            weights.insert(std::make_pair(weight, *iter));

            path.getLoopDetector().AddAlternative(iter->e_, weight);

        }
        NotifyAll(weights);

        EdgeContainer result;
        auto maxWeight = (--weights.end())->first;
        auto possibleEdge = weights.lower_bound(maxWeight / priorityCoefficient_);
        for (auto iter = possibleEdge; iter != weights.end(); ++iter) {
            if (wc_->IsExtensionPossible(iter->first)) {
                result.push_back(iter->second);
            }
        }
        return result;
    }
};



class MatePairExtensionChooser: public ExtensionChooser {

protected:

    void RemoveTrivial(BidirectionalPath& path) {
        wc_->GetExcludedEdges().clear();

        if (excludeTrivialWithBulges_) {
            analyzer_.ExcludeTrivialWithBulges(path, wc_->GetExcludedEdges());
        }
        else if (excludeTrivial_) {
            analyzer_.ExcludeTrivial(path, wc_->GetExcludedEdges());
        }
    }

public:
    MatePairExtensionChooser(Graph& g, WeightCounter * wc, double priority): ExtensionChooser(g, wc, priority) {

    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.empty()) {
            return edges;
        }

        RemoveTrivial(path);
        //TODO
        return edges;
    }
};



class ScaffoldingExtensionChooser: public ExtensionChooser {

    static bool compare(pair<int,double> a, pair<int,double> b)
    {
        if (a.first < b.first) return true;
        else return false;
    }

public:
	ScaffoldingExtensionChooser(Graph& g, WeightCounter * wc, double priority): ExtensionChooser(g, wc, priority) {

    }

	double AddInfoFromEdge(const std::vector<int>& distances, const std::vector<double>& weights, std::vector<pair<int,double> >& histogram, const BidirectionalPath& path, size_t j, double threshold)
	{
		double mean = 0.0;
		double sum  = 0.0;
		
		for (size_t l = 0; l < distances.size(); ++ l){
			if (distances[l] > max(0, (int) path.LengthAt(j) - (int) g_.k()) && weights[l] >= threshold) {
                mean += ((distances[l] - (int) path.LengthAt(j)) * weights[l]);
                sum += weights[l];
                histogram.push_back(make_pair(distances[l] - path.LengthAt(j), weights[l]));
			}
		}
		return mean / sum;
	}


    int CountMean(vector< pair<int,double> >& histogram)
    {
        double dist = 0.0;
        double sum = 0;
        for (size_t i = 0; i < histogram.size(); ++ i) {
			 sum += histogram[i].second;
		}
        for (size_t i = 0; i < histogram.size(); ++ i) {
			 dist += (histogram[i].first * histogram[i].second / sum);
		}
        return round(dist);
    }

    int CountDev(vector< pair<int,double> >& histogram, int avg)
    {
        double dev = 0.0;
        double sum = 0;
        for (size_t i = 0; i < histogram.size(); ++ i) {
             sum += histogram[i].second;
        }
        for (size_t i = 0; i < histogram.size(); ++ i) {
             dev += (((double) (histogram[i].first - avg)) * ((double) (histogram[i].first - avg)) * ((double) histogram[i].second));
        }
        return round(sqrt(dev / sum));
    }

    vector< pair<int,double> > FilterHistogram(vector< pair<int,double> >& histogram, int start, int end, int threshold)
    {
        vector< pair<int,double> > res;

        for (size_t i = 0; i < histogram.size(); ++i){
            if (histogram[i].first >= start && histogram[i].first <= end && histogram[i].second >= threshold) {
                res.push_back(histogram[i]);
            }
        }

        return res;
    }

    double CountAvrgDists(BidirectionalPath& path, EdgeId e, std::vector<pair<int,double> > & histogram)
    {
		std::vector<int> distances;
		std::vector<double> weights;
		
		double max_weight = 0.0;
		//bool print = true;
		for (size_t j = 0; j < path.Size(); ++ j) {
			wc_->GetDistances(path.At(j), e, distances, weights);
			
			for (size_t l = 0; l < weights.size(); ++ l){
				if (weights[l] > max_weight) {
				    max_weight = weights[l];
				}
			}

			if (distances.size() > 0) {
				AddInfoFromEdge(distances, weights, histogram, path, j, 0);

//				if (histogram.size() > 0) {
//					if (print) {
//						out << "\n************************" << endl;
//						out << "New pair :" << path.GetId() << " (" << g_.length(path[j])<< ")" << " - " <<  g_.int_id(e) << " (" << g_.length(e)<< ")" << endl;
//					}
//					out << "Mean = " << mean << endl;
//					print = false;
//				}
			}
			distances.clear();
		}
		return max_weight;
    }

    void FindBestFittedEdges(BidirectionalPath& path, EdgeContainer& edges, EdgeContainer& result)
    {
//    	ofstream out("./scaffolder.log", ostream::app);
//		out << "\n#########################################" << endl;
//		out << "Another path :" << path.GetId() << endl;
//		out << "Candidates:" << edges.size() << endl;

		std::vector<pair<int,double> > histogram;
		for (size_t i = 0; i < edges.size(); ++i){
			histogram.clear();
			double max_w = CountAvrgDists(path, edges[i].e_, histogram);

			for (int j = 0; j < 2; ++j) {
                int mean = CountMean(histogram);
                int dev = CountDev(histogram, mean);
                double cutoff = min(max_w * params.param_set.scaffolder_options.rel_cutoff, (double) params.param_set.scaffolder_options.cutoff);
//                if (mean != 0)
//                    out << "Mean: " << mean << " " << dev << endl;
                histogram = FilterHistogram(histogram, mean - (5 - j)  * dev, mean + (5 - j) * dev, cutoff);
			}

			double sum = 0.0;
			for (size_t j = 0; j < histogram.size(); ++j) {
			    sum += histogram[j].second;
			}

			if (sum > params.param_set.scaffolder_options.sum_threshold) {
				sort(histogram.begin(), histogram.end(), compare);
				int gap = CountMean(histogram);

//				out << "SUM: " << sum << endl;
//                out << "Histogram size: " << histogram.size() << endl;
//                out << "Avr gap value: " << gap << endl;
//                out << "weightI: " << wc_->CountIdealInfo(path, edges[i].e_, gap) << endl;
//                out << "weight1: " << wc_->CountWeight(path, edges[i].e_, gap) << endl;
//                out << "Histogram: \n";
//                for (size_t  ii = 0; ii < histogram.size(); ++ ii) {
//                    out << histogram[ii].first << " " << histogram[ii].second << endl;
//                }

                if (wc_->CountIdealInfo(path, edges[i].e_, gap) > 0.0) {
					result.push_back(EdgeWithDistance(edges[i].e_, gap));
					//out << "Added" << endl;
                }
			}
		}
//		out.close();
    }


    void FindBestFittedEdgesForClustered(BidirectionalPath& path, EdgeContainer& edges, EdgeContainer& result)
    {
        std::vector<pair<int,double> > histogram;
        for (size_t i = 0; i < edges.size(); ++i){
            histogram.clear();
            CountAvrgDists(path, edges[i].e_, histogram);

            double sum = 0.0;
            for (size_t j = 0; j < histogram.size(); ++j) {
                sum += histogram[j].second;
            }

            if (sum > params.param_set.scaffolder_options.cl_threshold) {
                sort(histogram.begin(), histogram.end(), compare);
                int gap = CountMean(histogram);

                if (wc_->CountIdealInfo(path, edges[i].e_, gap) > 0.0) {
                    result.push_back(EdgeWithDistance(edges[i].e_, gap));
                }
            }
        }
    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.empty()) {
            return edges;
        }
        EdgeContainer result;

        if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info) {
            FindBestFittedEdgesForClustered(path, edges, result);
        } else {
            FindBestFittedEdges(path, edges, result);
        }

        return result;
    }
};


}


#endif /* EXTENSION_HPP_ */
