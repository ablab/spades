/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file       BaseHypothesisEvaluator.cpp
//! @ingroup    SpadesHelpers
//! @brief      Combines code from the TS Basecaller and the TS Variant Caller to
//! @brief      give an indication about the feasibility of an alternative read sequence

#ifndef BASEHYPOTHESISEVALUATOR_H
#define BASEHYPOTHESISEVALUATOR_H

//#include "SpliceVariantHypotheses.h"

//#include "BaseCallerUtils.h"
#include "TreephaserLite.h"
#include <bamtools/api/BamAlignment.h>

using namespace std;
using namespace ion;
using namespace BamTools;


// Function to calculate signal predictions
//! @brief  Function to calculate signal predictions and evaluate L2-norm.
//! @param[in]  alignment      A read in the BamAlignment format
//! @param[in]  flow_order_str Flow order string
//! @param[in]  alt_base_hyp   Alternative base hypothesis
//! @param[out] delta_score    Difference in fit between hypotheses; negative number indicate better alt.
//! @param[out] fit_score      Overall fit, a higher number indicates a noisier read or a problem.
void BaseHypothesisEvaluator(BamTools::BamAlignment    &alignment,
                             const string              &flow_order_str,
                             const string              &alt_base_hyp,
                             float                     &delta_score,
                             float                     &fit_score,
                             int                       heavy_verbose=0);


// Retrieve basecaller specific bam tags from BamTools alignment
bool GetBamTags(BamTools::BamAlignment &alignment,
		                const int              &num_flows,
                        vector<float>          &measurements,
		                vector<float>          &phase_params,
                        int                    &start_flow);


// Solve for hard and soft clipped bases at the start of the read, before start_flow
int GetMasterReadPrefix(TreephaserLite       &treephaser,
		                const ion::FlowOrder &flow_order,
                        const int            &start_flow,
                        const string         &called_bases,
                        BasecallerRead       &master_read);


// Print out some messages
void PredictionGenerationVerbose(const vector<string>         &Hypotheses,
                                 const vector<BasecallerRead> &hypothesesReads,
                                 const vector<float>          &phase_params,
                                 const ion::FlowOrder         &flow_order,
                                 const int                    &start_flow,
                                 const int                    &prefix_size);

// Complement a nucleotide
char NucComplement(char nuc);

// Reverse complement a base string in place
void RevComplementInPlace(string& seq);

#endif // BASEHYPOTHESISEVALUATOR_H
