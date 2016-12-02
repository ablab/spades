//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../online_vis/environment.hpp"
#include "compare_standard.hpp"
#include "coloring.hpp"
#include "coordinates_handler.hpp"
#include "test_utils.hpp"
#include "serialization.hpp"

namespace online_visualization {

class CapEnvironmentManager;
/*
 * Cap Environment is designed to handle operations on fixed
 * set of genomes
 */
class CapEnvironment : public Environment {
  friend class CapEnvironmentManager;

 private:
  typedef debruijn_graph::Graph Graph;
  typedef Graph::VertexId VertexId;
  typedef Graph::EdgeId EdgeId;

  typedef debruijn_graph::KmerStoringEdgeIndex<Graph, RtSeq, kmer_index_traits<RtSeq>, debruijn_graph::SimpleStoring> RtSetIndex;
  typedef debruijn_graph::graph_pack<Graph, RtSeq, RtSetIndex> RtSeqGraphPack;
  typedef debruijn_graph::KmerStoringEdgeIndex<Graph, cap::LSeq, kmer_index_traits<cap::LSeq>, debruijn_graph::SimpleStoring> LSeqIndex;
  typedef debruijn_graph::graph_pack<Graph, cap::LSeq, LSeqIndex> LSeqGraphPack;

  typedef cap::ColorHandler<Graph> ColorHandler;
  typedef cap::CoordinatesHandler<Graph> CoordinatesHandler;

  std::string name_;
  std::string dir_;
  std::string description_;

  // Sequence of Ks which were used to refine current genome
  std::vector<unsigned> k_history_;
  // History of number of studied genomes (for case of adding genomes in the middle of pipeline)
  std::vector<size_t> num_genomes_history_;

  // Paths on fs
  std::vector<std::string> init_genomes_paths_;
  // Genome sequences themselves. Yes, it may be lots of GBs.
  std::vector<Sequence> genomes_;
  std::vector<std::string> genomes_names_;

  std::shared_ptr<RtSeqGraphPack> gp_rtseq_;
  std::shared_ptr<LSeqGraphPack> gp_lseq_;

  std::shared_ptr<ColorHandler> coloring_;

  CoordinatesHandler coordinates_handler_;

  // Aliases to GraphPack parts:
  //
  Graph *graph_;
  EdgesPositionHandler<Graph> *edge_pos_;
  GraphElementFinder<Graph> *element_finder_;

  // Information concerning the default way to write out info about diversities
  std::string event_log_path_;
  // Either "w" or "a", and using the "a" mode during the environment load file
  // will be anyway recreated (and purged) though.
  std::string event_log_file_mode_;

  // Environment Manager for complex methods on this Environment
  std::shared_ptr<CapEnvironmentManager> manager_;

  void AssignGPReferences() {
    VERIFY(gp_lseq_ != NULL || gp_rtseq_ != NULL);
    if (LSeqIsUsed()) {
      graph_ = &(gp_lseq_->g);
      edge_pos_ = &(gp_lseq_->edge_pos);
      element_finder_ = &(gp_lseq_->element_finder);
    } else {
      graph_ = &(gp_rtseq_->g);
      edge_pos_ = &(gp_rtseq_->edge_pos);
      element_finder_ = &(gp_rtseq_->element_finder);
    }
    coordinates_handler_.SetGraph(graph_);
  }

  void set_gp(const std::shared_ptr<LSeqGraphPack> &gp_lseq) {
    gp_lseq_ = gp_lseq;
  }

  void set_gp(const std::shared_ptr<RtSeqGraphPack> &gp_rtseq) {
    gp_rtseq_ = gp_rtseq;
  }

 public:
  static const unsigned kNoGraphK = -1;
  const std::string kDefaultGPWorkdir;

  CapEnvironment(const std::string &name, /*const std::string base_path, */const std::string &description = "")
      : Environment(name, cap_cfg::get().cache_root + "/env_" + name),
        name_(name),
        dir_(cap_cfg::get().cache_root + "/env_" + name),
        //base_path_(base_path),
        description_(description),
        k_history_(),
        num_genomes_history_(),
        init_genomes_paths_(),
        gp_rtseq_(),
        gp_lseq_(),
        coloring_(),
        coordinates_handler_(),
        graph_(NULL),
        edge_pos_(NULL),
        element_finder_(NULL),
        event_log_path_(dir_ + "/" + cap_cfg::get().default_log_filename),
        event_log_file_mode_(cap_cfg::get().default_log_file_mode),
        manager_(std::make_shared<CapEnvironmentManager>(this)),
        kDefaultGPWorkdir("./tmp") {
    cap::utils::MakeDirPath(dir_);
  }

  void WriteToStream(std::ostream &out) const {
    cap::Serializer s(out);

    s.WriteLine("name", name_);
    s.WriteLine("dir", dir_);
    s.WriteLine("description", description_);
    s.WriteLine("k_history", k_history_);
    s.WriteLine("num_genomes_history", num_genomes_history_);
    s.WriteLine("init_genomes_paths", init_genomes_paths_);
    s.WriteLine("genomes_names", genomes_names_);

    s.WriteLine("genomes", genomes_);

    s.WriteLine("coordinates_threads", coordinates_handler_.GetStoredThreads());
  }

  void ReadFromStream(std::istream &in) {
    cap::Deserializer s(in);

    s.ReadStream();

    s.ReadValue("name", name_);
    s.ReadValue("dir", dir_);
    s.ReadValue("description", description_);
    s.ReadValue("k_history", k_history_);
    s.ReadValue("num_genomes_history", num_genomes_history_);
    s.ReadValue("init_genomes_paths", init_genomes_paths_);
    s.ReadValue("genomes_names", genomes_names_);

    s.ReadValue("genomes", genomes_);

    std::vector<std::pair<uint, std::vector<CoordinatesHandler::Thread>>> coords_threads;
    s.ReadValue("coordinates_threads", coords_threads);
    coordinates_handler_.SetStoredThreads(coords_threads);
  }

  void CheckConsistency() const {
    VERIFY(gp_rtseq_ == NULL || gp_lseq_ == NULL);
    bool have_any_gp = gp_rtseq_ != NULL || gp_lseq_ != NULL;
    if (have_any_gp) {
      VERIFY(graph_ != NULL);
      VERIFY(edge_pos_ != NULL);
      VERIFY(element_finder_ != NULL);
      VERIFY(coordinates_handler_.GetGraph() == graph_);
    }
  }

  void ClearGP() {
    // shared_ptr deletes automatically
    coordinates_handler_.UnsetGraph();
    gp_rtseq_.reset();
    gp_lseq_.reset();
    graph_ = NULL;
    edge_pos_ = NULL;
    element_finder_ = NULL;
    coloring_.reset();

    CheckConsistency();
  }

  template <class GraphPack>
  void SetGraphPack(const std::shared_ptr<GraphPack> &gp) {
    VERIFY(gp_rtseq_ == NULL && gp_lseq_ == NULL);
    set_gp(gp);
    AssignGPReferences();
    CheckConsistency();
  }

  unsigned GetGraphK() const {
    if (gp_rtseq_ == NULL && gp_lseq_ == NULL) {
      return kNoGraphK;
    }

    return unsigned(graph_->k());
  }
  bool LSeqIsUsed() const {
    return gp_lseq_ != NULL;
  }

  // Method defining of we use RtSeq or LSeq for particular K
  bool UseLSeqForThisK(unsigned k) const {
    return k > 201;
  }

  const std::string &name() const {
    return name_;
  }
  const std::string &dir() const {
    return dir_;
  }
  const Graph &graph() const {
    return *graph_;
  }

  RtSeqGraphPack& rt_seq_gp() const {
      return *gp_rtseq_;
  }

  LSeqGraphPack& l_seq_gp() const {
      return *gp_lseq_;
  }

  const vector<Sequence>& genomes() const {
      return genomes_;
  }

  vector<Sequence>& mutable_genomes() {
      return genomes_;
  }

  const vector<string>& genome_names() const {
      return genomes_names_;
  }

  const CoordinatesHandler& coordinates_handler() const {
      return coordinates_handler_;
  }

  unsigned genome_cnt() const {
      return unsigned(genomes_names_.size());
  }

  const EdgesPositionHandler<Graph> &edge_pos() const {
    return *edge_pos_;
  }
  const ColorHandler &coloring() const {
    return *coloring_;
  }
  const std::string &event_log_path() {
    return event_log_path_;
  }
  const std::string &event_log_file_mode() {
    return event_log_file_mode_;
  }
  CapEnvironmentManager &manager() const {
    return *manager_;
  }

};
}
