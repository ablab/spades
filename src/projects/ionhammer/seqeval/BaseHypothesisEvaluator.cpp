/* Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved */

//! @file       BaseHypothesisEvaluator.cpp
//! @ingroup    SpadesHelpers
//! @brief      Combines code from the TS Basecaller and the TS Variant Caller to
//! @brief      give an indication about the feasibility of an alternative read sequence

#include "BaseHypothesisEvaluator.h"

// Function to fill in predicted signal values
void BaseHypothesisEvaluator(BamTools::BamAlignment    &alignment,
                             const string              &flow_order_str,
                             const string              &alt_base_hyp,
                             float                     &delta_score,
                             float                     &fit_score,
                             int                       heavy_verbose) {

    // --- Step 1: Initialize Objects and retrieve relevant tags

    delta_score = 1e5;
    fit_score   = 1e5;
    vector<string>   Hypotheses(2);
    vector<float>    measurements, phase_params;
    int              start_flow, num_flows, prefix_flow=0;

    if (not GetBamTags(alignment, flow_order_str.length(), measurements, phase_params, start_flow))
      return;
    num_flows = measurements.size();
    ion::FlowOrder flow_order(flow_order_str, num_flows);
    BasecallerRead master_read;
    master_read.SetData(measurements, flow_order.num_flows());
    TreephaserLite   treephaser(flow_order);
    treephaser.SetModelParameters(phase_params[0], phase_params[1]);

    // --- Step 2: Solve beginning of the read
    // Look at mapped vs. unmapped reads in BAM
    Hypotheses[0] = alignment.QueryBases;
    Hypotheses[1] = alt_base_hyp;
    // Safety: reverse complement reverse strand reads in mapped bam
    if (alignment.IsMapped() and alignment.IsReverseStrand()) {
      RevComplementInPlace(Hypotheses[0]);
      RevComplementInPlace(Hypotheses[1]);
    }

    prefix_flow = GetMasterReadPrefix(treephaser, flow_order, start_flow, Hypotheses[0], master_read);
    unsigned int prefix_size = master_read.sequence.size();

    // --- Step 3: creating predictions for the individual hypotheses

    vector<BasecallerRead> hypothesesReads(Hypotheses.size());
    vector<float> squared_distances(Hypotheses.size(), 0.0);
    int max_last_flow = 0;

    for (unsigned int i_hyp=0; i_hyp<hypothesesReads.size(); ++i_hyp) {

      hypothesesReads[i_hyp] = master_read;
      // --- add hypothesis sequence to clipped prefix
      unsigned int i_base = 0;
      int i_flow = prefix_flow;

      while (i_base<Hypotheses[i_hyp].length() and i_base<(2*(unsigned int)flow_order.num_flows()-prefix_size)) {
        while (i_flow < flow_order.num_flows() and flow_order.nuc_at(i_flow) != Hypotheses[i_hyp][i_base])
          i_flow++;
        if (i_flow < flow_order.num_flows() and i_flow > max_last_flow)
          max_last_flow = i_flow;
        if (i_flow >= flow_order.num_flows())
          break;
        // Add base to sequence only if it fits into flow order
        hypothesesReads[i_hyp].sequence.push_back(Hypotheses[i_hyp][i_base]);
        i_base++;
      }
      i_flow = min(i_flow, flow_order.num_flows()-1);

      // Solver simulates beginning of the read and then fills in the remaining clipped bases for which we have flow information
      treephaser.Solve(hypothesesReads[i_hyp], num_flows, i_flow);
    }
    // Compute L2-distance of measurements and predictions
    for (unsigned int i_hyp=0; i_hyp<hypothesesReads.size(); ++i_hyp) {
      for (int iFlow=0; iFlow<=max_last_flow; iFlow++)
        squared_distances[i_hyp] += (measurements.at(iFlow) - hypothesesReads[i_hyp].prediction.at(iFlow)) *
                                    (measurements.at(iFlow) - hypothesesReads[i_hyp].prediction.at(iFlow));
    }

    // Delta: L2-distance of alternative base Hypothesis - L2-distance of bases as called
    delta_score = squared_distances.at(1) - squared_distances.at(0);
    fit_score   = min(squared_distances.at(1), squared_distances.at(0));


    // --- verbose ---
    if (heavy_verbose > 1 or (delta_score < 0 and heavy_verbose > 0)) {
      cout << "Processed read " << alignment.Name << endl;
      cout << "Delta Fit: " << delta_score << " Overall Fit: " << fit_score << endl;
      PredictionGenerationVerbose(Hypotheses, hypothesesReads, phase_params, flow_order, start_flow, prefix_size);
    }

}

// ----------------------------------------------------------------------

bool GetBamTags(BamTools::BamAlignment &alignment,
                        const int              &num_flows,
                        vector<float>          &measurements,
                        vector<float>          &phase_params,
                        int                    &start_flow) {

  vector<int16_t>  quantized_measurements;
  // Retrieve normalized measurements from BAM file
  if (not alignment.GetTag("ZM", quantized_measurements)) {
    cerr << "ERROR: Normalized measurements ZM:tag is not present in read " << alignment.Name << endl;
    return false;
  }
  if ((int)quantized_measurements.size() > num_flows) {
    cerr << "ERROR: Normalized measurements ZM:tag length exceeds flow order length in read " << alignment.Name << endl;
    return false;
  }
  measurements.assign(quantized_measurements.size(), 0.0);
  for (size_t counter = 0; counter < quantized_measurements.size(); ++counter)
    measurements.at(counter) = (float)quantized_measurements.at(counter)/256;

  // Retrieve phasing parameters from BAM file
  if (not alignment.GetTag("ZP", phase_params)) {
    cerr << "ERROR: Phasing Parameters ZP:tag is not present in read " << alignment.Name << endl;
    return false;
  }
  if (phase_params.size() != 3) {
    cerr << "ERROR: Phasing Parameters ZP:tag does not have 3 phase parameters in read " << alignment.Name << endl;
    return false;
  }
  if (phase_params[0] < 0 or phase_params[0] > 1 or phase_params[1] < 0 or phase_params[1] > 1
      or phase_params[2] < 0 or phase_params[2] > 1) {
    cerr << "ERROR: Phasing Parameters ZP:tag outside of [0,1] range in read " << alignment.Name << endl;
    return false;
  }
  phase_params[2] = 0.0f;   // ad-hoc corrector: zero droop

  // Retrieve start flow
  if (not alignment.GetTag("ZF", start_flow)) {
    cerr << "ERROR: Start Flow ZF:tag not found in read " << alignment.Name << endl;
    return false;
  }
  if (start_flow < 0 or start_flow >= num_flows) {
    cerr << "ERROR: Start flow outsize of [0,num_flows) range in read " << alignment.Name << endl;
    cerr << "Start flow: " << start_flow << " Number of flows: " << num_flows;
    return false;
  }
  // A start flow of zero indicated a read that did not pass basecaller filters
  if (start_flow == 0) {
    cerr << "WARNING: Start Flow ZF:tag has zero value in read " << alignment.Name << endl;
    return false;
  }
  return true;
}

// ----------------------------------------------------------------------

int GetMasterReadPrefix(TreephaserLite       &treephaser,
                        const ion::FlowOrder &flow_order,
                        const int            &start_flow,
                        const string         &called_bases,
                        BasecallerRead       &master_read) {

  // Solve beginning of maybe clipped read
  int until_flow = min((start_flow+20), flow_order.num_flows());
  treephaser.Solve(master_read, until_flow, 0);

  // StartFlow clipped? Get solved HP length at startFlow.
  unsigned int base = 0;
  int flow = 0;
  unsigned int HPlength = 0;
  while (base < master_read.sequence.size()) {
    while (flow < flow_order.num_flows() and flow_order.nuc_at(flow) != master_read.sequence[base]) {
      flow++;
    }
    if (flow > start_flow or flow == flow_order.num_flows())
      break;
    if (flow == start_flow)
      HPlength++;
    base++;
  }
  //if (global_context.DEBUG>2)
  //  printf("Solved %d bases until (not incl.) flow %d. HP of height %d at flow %d.\n", base, flow, HPlength, start_flow);

  // Get HP size at the start of the read as called in Hypotheses[0]
  unsigned int count = 1;
  while (count < called_bases.length() and called_bases.at(count) == called_bases.at(0))
    count++;
  //if (global_context.DEBUG>2)
  //  printf("Hypothesis starts with an HP of length %d\n", count);
  // Adjust the length of the prefix and erase extra solved bases
  if (HPlength>count)
    base -= count;
  else
    base -= HPlength;
  master_read.sequence.erase(master_read.sequence.begin()+base, master_read.sequence.end());

  // Get flow of last prefix base
  int prefix_flow = 0;
  for (unsigned int i_base = 0; i_base < master_read.sequence.size(); i_base++) {
    while (prefix_flow < flow_order.num_flows() and flow_order.nuc_at(prefix_flow) != master_read.sequence[i_base])
      prefix_flow++;
  }

  return prefix_flow;
}


// ----------------------------------------------------------------------

void PredictionGenerationVerbose(const vector<string>         &Hypotheses,
                                 const vector<BasecallerRead> &hypothesesReads,
                                 const vector<float>          &phase_params,
                                 const ion::FlowOrder         &flow_order,
                                 const int                    &start_flow,
                                 const int                    &prefix_size) {

  printf("Calculating predictions for %d hypotheses starting at flow %d:\n", (int)Hypotheses.size(), start_flow);
  for (unsigned int iHyp=0; iHyp<Hypotheses.size(); ++iHyp) {
    for (unsigned int iBase=0; iBase<Hypotheses[iHyp].length(); ++iBase)
      printf("%c", Hypotheses[iHyp][iBase]);
    printf("\n");
  }
  printf("Solved read prefix: ");
  for (int iBase=0; iBase<prefix_size; ++iBase)
    printf("%c", hypothesesReads[0].sequence[iBase]);
  printf("\n");
  printf("Extended Hypotheses reads to:\n");
  for (unsigned int iHyp=0; iHyp<hypothesesReads.size(); ++iHyp) {
    for (unsigned int iBase=0; iBase<hypothesesReads[iHyp].sequence.size(); ++iBase)
      printf("%c", hypothesesReads[iHyp].sequence[iBase]);
    printf("\n");
  }
  printf("Phasing Parameters, cf: %f ie: %f dr: %f \n Predictions: \n",
          phase_params[0], phase_params[1], phase_params[2]);
  cout << "Flow Order  : ";
  for (int i_flow=0; i_flow<flow_order.num_flows(); i_flow++) {
    cout << flow_order.nuc_at(i_flow) << "    ";
    if (hypothesesReads[0].normalized_measurements[i_flow] < 0)
      cout << " ";
  }
  cout << endl << "Flow Index  : ";
  for (int i_flow=0; i_flow<flow_order.num_flows(); i_flow++) {
      cout << i_flow << " ";
      if (i_flow<10)        cout << "   ";
      else if (i_flow<100)  cout << "  ";
      else if (i_flow<1000) cout << " ";
      if (hypothesesReads[0].normalized_measurements[i_flow] < 0)
        cout << " ";
    }
  cout << endl << "Measured    : ";
  for (unsigned int i_flow=0; i_flow<hypothesesReads[0].normalized_measurements.size(); ++i_flow) {
    printf("%.2f", hypothesesReads[0].normalized_measurements[i_flow]);
    if (hypothesesReads[0].normalized_measurements[i_flow] < 10)
      cout << " ";
  }
  cout << endl;
  for (unsigned int i_Hyp=0; i_Hyp<hypothesesReads.size(); ++i_Hyp) {
    cout << "Prediction "<< i_Hyp << ": ";
    for (unsigned int i_flow=0; i_flow<hypothesesReads[i_Hyp].prediction.size(); ++i_flow) {
      printf("%.2f", hypothesesReads[i_Hyp].prediction[i_flow]);
      if (hypothesesReads[i_Hyp].prediction[i_flow] < 10)
        cout << " ";
      if (hypothesesReads[0].normalized_measurements[i_flow] < 0)
        cout << " ";
    }
    cout << endl;
  }
  cout << " ------------------- " << endl;
}

// ----------------------------------------------------------------------

char NucComplement (char nuc)
{
  switch(nuc) {
    case ('A') : return 'T';
    case ('C') : return 'G';
    case ('G') : return 'C';
    case ('T') : return 'A';
    case ('a') : return 't';
    case ('c') : return 'g';
    case ('g') : return 'c';
    case ('t') : return 'a';

    default:  return nuc; // e.g. 'N' and '-' handled by default
  }
}

void RevComplementInPlace(string& seq) {

  char c;
  int forward_idx = 0;
  int backward_idx = seq.size()-1;
  while (forward_idx < backward_idx) {
    c = seq[forward_idx];
    seq[forward_idx]  = NucComplement(seq[backward_idx]);
    seq[backward_idx] = NucComplement(c);
    forward_idx++;
    backward_idx--;
  }
  if (forward_idx == backward_idx)
    seq[forward_idx] = NucComplement(seq[forward_idx]);
}
