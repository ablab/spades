/**
 * bayes_quality.h
 *
 *  Created on: Mar 23, 2011
 *      Author: snikolenko
 */
#ifndef BAYES_QUALITY_H_
#define BAYES_QUALITY_H_

#include <vector>
#include "sequence.hpp"
#include "logging.hpp"

using namespace std;

namespace bayes_quality {


// Vector of quality values
typedef vector<int> QVector;
// A read with quality values
typedef pair<Sequence, QVector> QRead;

/**
 * @brief Vertex of a condensed graph.
 *
 * Class representing vertex of condensed graph.
 * Contains the Sequence of nucleotides, average coverage of this sequence,
 * pointers to children and coverage of corresponding edges.
 *
 */
class BayesQualityGenome {
private:
	Sequence genome_;

public:
	//bool deleted;
	BayesQualityGenome(const char *genome) :
		genome_(genome) { }

	double ReadBQInt(QRead);
};

}

#endif /* BAYES_QUALITY_H_ */

