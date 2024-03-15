/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file     TreephaserLite.h
//! @ingroup  Spades Helper
//! @brief    A lighter version of TS-TreephaserLite. Perform dephasing and call base sequence by tree search.

#ifndef TREEPHASERLITE_H
#define TREEPHASERLITE_H

// Maximum length of a Homopolymer, signal processing does not go any higher
#define MAX_HPXLEN 23

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <stddef.h>
#include <algorithm>
#include "BaseCallerUtils.h"

using namespace std;


//! @brief    Input/output data structure for TreephaserLite
//! @ingroup  BaseCaller

struct BasecallerRead {

  void SetData(const vector<float> &measurements, int num_flows);
  void SetDataAndKeyNormalize(const float *measurements, int num_flows, const int *key_flows, int num_key_flows);

  float           key_normalizer;           //!< Scaling factor used for initial key normalization
  vector<float>   raw_measurements;         //!< Measured, key-normalized flow signal
  vector<float>   normalized_measurements;  //!< Measured flow signal with best normalization so far
  vector<float>   prediction;               //!< Model-based phased signal predicted for the "solved" sequence
  vector<char>    sequence;                 //!< Vector of ACGT bases. Output from Solver, input to Simulator
};


const  int    kMinWindowSize_     = 20;   //!< Minimum normalization window size
const  int    kMaxWindowSize_     = 60;   //!< Maximum normalization window size


//! @brief    Performs dephasing and base calling by tree search
//! @ingroup  BaseCaller
//! @details
//! TreephaserLite is responsible for determining base sequence from phased signal output by background model.
//! It uses a generative phasing model that can produce expected signal (prediction) for a partial
//! or complete base sequence. It further uses tree search to find a base sequence with prediction
//! matching the background model signal (measurements) the closest.
//! Additionally, TreephaserLite contains signal normalization procedures that can correct for additive
//! and multiplicative distortion using earlier predicted signal.
//! This allows dephasing and normalization to be performed iteratively (adaptive normalization)

class TreephaserLite {

public:
    // These need to be public for TreephaserSSE to use.
    const static int    kWindowSizeDefault_ = 38;   //!< Default normalization window size

  //! @brief  Constructor.
  //! @param[in] flow_order   Flow order object, also stores number of flows
  TreephaserLite(const ion::FlowOrder& flow_order, const int windowSize=kWindowSizeDefault_);

  //! @brief  Set the normalization window size
  //! @param[in]  windowSize  Size of the normalization window to use.
  inline void SetNormalizationWindowSize(const int windowSize)
      { windowSize_ = max(kMinWindowSize_, min(windowSize, kMaxWindowSize_));}

  //! @brief  Initializes phasing model using specific phasing parameters.
  //! @param[in]  cf          Carry forward rate, how much nuc from previous flow is encountered
  //! @param[in]  ie          Incomplete extension rate, how much polymerase fails to incorporate
  //! @param[in]  dr          Droop, how much polymerase deactivates during an incorporation
  void  SetModelParameters(double cf, double ie, double dr=0.0);

  //! @brief  Perform adaptive normalization using WindowedNormalize
  //! @param[in,out]  read           Input and output information for the read
  //! @param[in]      max_flows      Number of flows to process
  //! @param[in]      sliding_window Switch to use a sliding window in normalization
  void  NormalizeAndSolve(BasecallerRead& read, int max_flows, bool sliding_window=true);

  //! @brief  Tree-search-based dephasing.
  //! @param[in]  read.normalized_measurements    Normalized measurements
  //! @param[out] read.sequence   Predicted base sequence
  //! @param[out] read.prediction Predicted signal
  //! @param[in]  max_flows       Number of flows to process
  //! @param[in]  restart_flows   Number of flows to simulate, rather than solve
  void  Solve(BasecallerRead& read, int max_flows, int restart_flows=0);

  //! @brief  Generate predicted signal from base sequence
  //! @param[in]  read.sequence     Base sequence
  //! @param[out] read.prediction   Predicted signal
  //! @param[in]  max_flows         Number of flows to process
  void  Simulate(BasecallerRead& read, int max_flows);

  //! @brief  Correct for flow-varying multiplicative and additive distortion
  //! @param[in]  read.prediction               Model-predicted signal
  //! @param[in]  read.raw_measurements         Flow signal before normalization
  //! @param[out] read.normalized_measurements  Flow signal after normalization
  //! @param[in]  num_steps                     Number of windows-worth of predictions to use
  //! @param[in]  window_size                   Size of a window in flows
  void  WindowedNormalize(BasecallerRead& read, int num_steps, int window_size) const;


  //! @brief    Treephaser's slot for partial base sequence, complete with tree search metrics and state for extending
  struct TreephaserPath {
    bool              in_use;                   //!< Is this slot in use?

    // Phasing and input-output state of this path
    int               flow;                     //!< In phase flow of last incorporated base
    vector<float>     state;                    //!< Histogram of flows at which last base was incorporated
    int               window_start;             //!< Start flow (inclusive) of meaningful state values
    int               window_end;               //!< End flow (noninclusive) of meaningful state values
    vector<float>     prediction;               //!< Model-based phased signal predicted for this path
    vector<char>      sequence;                 //!< Vector of ACGT bases corresponding to this path
    int               last_hp;                  //!< Length of the last homopolymer in sequence

    // Path metrics and related values
    float             path_metric;              //!< Primary tree search metric
    float             residual_left_of_window;  //!< Residual left of the state window
    float             per_flow_metric;          //!< Auxiliary tree search metric, useful for stack pruning
    int               dot_counter;              //!< Number of extreme mismatch flows encountered so far
  };

  TreephaserPath& path(int idx) { return path_[idx]; }


  //! @brief  Set path to an empty sequence, a starting point for phasing simulation
  //! @param[out]  state    Path slot
  void InitializeState(TreephaserPath *state) const;

  //! @brief  Perform a path extension by one nucleotide
  //! @param[out]  child     Path slot to store the extended path
  //! @param[in]   parent    Path to be extended
  //! @param[in]   nuc       Nucleotide (integer) to extend the path by
  //! @param[in]   max_flow  Do not read/write past this flow
  void AdvanceState(TreephaserPath *child, const TreephaserPath *parent, char nuc, int max_flow) const;

  //! @brief  Perform a path extension by one nucleotide
  //! @param[in,out] state     Path to be extended in place
  //! @param[in]     nuc       Nucleotide (integer) to extend the path by
  //! @param[in]     max_flow  Do not read/write past this flow
  void AdvanceStateInPlace(TreephaserPath *state, char nuc, int max_flow) const;


protected:

  int                 windowSize_;                //!< Normalization window size

  ion::FlowOrder      flow_order_;                //!< Sequence of nucleotide flows
  vector<float>       transition_base_[8];        //!< Probability of polymerase incorporating and staying active
  vector<float>       transition_flow_[8];        //!< Probability of polymerase not incorporating and staying active
  vector<TreephaserPath> path_;                   //!< Preallocated space for partial path slots

  // Magic constants
  constexpr static int    kNumPaths = 8;              //!< Maximum number of paths considered
  constexpr static float  kExtendThreshold = 0.2f;     //!< Threshold for extending paths
  constexpr static float  kNegativeMultiplier = 2.0f;  //!< Extra weight on the negative residuals
  constexpr static float  kDotThreshold = 0.3f;        //!< percentage of expected Signal that constitutes a "dot"
  constexpr static int    kMaxHP = MAX_HPXLEN;        //!< Maximum callable homopolymer length
  constexpr static float  kStateWindowCutoff = 1e-6f;  //!< Minimum fraction to be taken into account
  constexpr static int    kMaxPathDelay = 40;         //!< Paths that are delayed more are killed
};

#endif // TREEPHASERLITE_H
