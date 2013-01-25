//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * analyzer.hpp
 *
 *  Created on: Mar 29, 2012
 *      Author: andrey
 */

#ifndef ANALYZER_HPP_
#define ANALYZER_HPP_

#include "pe_config_struct.hpp"
#include "extension_chooser.hpp"

namespace path_extend {

class ResolveAnalyzer: public ExtensionChooserListener {

protected:
    std::vector <double> pathWeights;

    std::vector < std::vector<double> > pathAlternatives;

public:

    virtual void ExtensionChosen(double weight) {
        pathWeights.push_back(weight);
    }

    virtual void ExtensionChosen(AlternativeConteiner& alts) {
        std::vector<double> weights;

        for (auto iter = alts.rbegin(); iter != alts.rend(); ++iter) {
            weights.push_back(iter->first);
        }
        pathAlternatives.push_back(weights);
        pathWeights.push_back(weights.back());
    }

    void PrintWeightHisogram(double h = 0.1) {
        VERIFY_MSG(math::eq(h, 0.0), "Histogram step cannot be zero")

        double maxWeight = 0.0;
        for (size_t i = 0; i < pathWeights.size(); ++i) {
            if (math::gr(pathWeights[i], maxWeight)) {
                maxWeight = pathWeights[i];
            }
        }

        std::vector<int> hist;
        hist.resize(((int) math::round(maxWeight / h)) + 1);
        for (size_t i = 0; i < pathWeights.size(); ++i) {
            int pos = math::round(pathWeights[i] / h);
            if (pos >= 0 && pos < hist.size()) {
                ++hist[pos];
            }
        }

        INFO("Weight histogram");
        for (size_t i = 0; i < hist.size(); ++i) {
            INFO(h * i << " " << hist[i]);
        }
    }

};

}


#endif /* ANALYZER_HPP_ */
