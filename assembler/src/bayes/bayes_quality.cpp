/*
 * bayes_quality.cpp
 *
 *  Created on: 23.03.2011
 *      Author: snikolenko
 */

#include "bayes_quality.hpp"
#include <algorithm>
#include <numeric>
#include <math.h>

using namespace std;

LOGGER("a");

namespace bayes_quality {

double BayesQualityGenome::ReadBQInt(QRead r) {

	// initialize temporary vector of possible start positions
	QVector qv(genome_.size() - r.first.size());
	
	// compute the sum of all quality values
	int totalQ = accumulate(r.second.begin(), r.second.end(), 0);
	
	// fill the vector with that sum
	fill(qv.begin(), qv.end(), totalQ);
	
	// for each element of the read
	for (size_t i = 0; i < r.first.size(); ++i) {
		// subtract its quality value if the corresponding elements match
		for (size_t j = 0; j < qv.size(); ++j) {
			if ( r.first[i] == genome_[i+j] ) {
				qv[j] -= r.second[i];
			}
		}
	}
	
	// now add them all up
	double res = 0;
	for (size_t i = 0; i < qv.size(); ++i) {
		res += pow(10, -qv[i]/10.0);
	}
	
	// return the sum of the qv vector
	return res;
}

}

