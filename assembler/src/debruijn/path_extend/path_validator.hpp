/*;
 * paths_validator.hpp
 *
 *  Created on: Jul 4, 2012
 *      Author: tanunia
 */

#ifndef PATH_VALIDATOR_HPP_
#define PATH_VALIDATOR_HPP_

#include "bidirectional_path.hpp"
#include <cstdlib>

using namespace path_extend;

class PathValidator{
	Graph& g_;
	WeightCounter* wc_;
	WeightCounter* swc_;
	double threshold_;

public:

    PathValidator(Graph& g, WeightCounter * wc, WeightCounter * swc, double threshold = 0.0): g_(g), wc_(wc), swc_(swc), threshold_(threshold){
    }

	bool ValidateSinglePath(BidirectionalPath& path)
	{
		bool isValid = true;
		double weight;
        for (size_t i = 0; i < path.Size(); ++i) {
        	EdgeId e = path[i];
        	if (path.GapAt(i) == 0) {
        		weight = wc_->CountWeight(path, e);
        	} else {
        		weight = swc_->CountWeight(path, e);
        	}
        	if (weight - threshold_ < 0) {
        		 isValid=false;
        		 break;
        	}
        }
		return isValid;
	}

	std::vector<bool> ValidatePaths(PathContainer& paths)
	{
		std::vector<bool> validationResult(paths.size());
		for (size_t i = 0; i < paths.size(); ++i) {
			validationResult[i] = ValidateSinglePath(*paths.Get(i));
		}
		return validationResult;
	}

	std::vector<bool> ValidatePaths(std::vector < std::vector < EdgeId > > & paths)
	{
		std::vector<bool> validationResult(paths.size());
		for (size_t i = 0; i < paths.size(); ++i) {
			BidirectionalPath path = BidirectionalPath(g_,paths[i]);
			validationResult[i] = ValidateSinglePath(path);
		}
		return validationResult;
	}

};

#endif /* PATH_VALIDATOR_HPP_ */



