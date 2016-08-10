//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "stages/simplification_pipeline/graph_simplification.hpp"

#include "compare_standard.hpp"
#include "pipeline/graphio.hpp"

#include "comparison_utils.hpp"
#include "diff_masking.hpp"
#include "genome_correction.hpp"
#include "assembly_compare.hpp"
#include "simple_indel_finder.hpp"
#include "simple_inversion_finder.hpp"
#include "graph_traversal_constraints.hpp"
#include "test_utils.hpp"

#include "cap_environment.hpp"
#include "io/reads/sequence_reader.hpp"
#include "pipeline/config_struct.hpp"
#include "junk_cropping_reader.hpp"

namespace online_visualization {

using namespace cap;

class CapEnvironmentManager {
  typedef CapEnvironment::Graph Graph;
  typedef CapEnvironment::LSeqGraphPack LSeqGP;
  typedef CapEnvironment::RtSeqGraphPack RtSeqGP;
  typedef Graph::VertexId VertexId;
  typedef Graph::EdgeId EdgeId;

  CapEnvironment *env_;
  //vector<ContigStream *> last_streams_used_;

  void WriteStateDesc(std::string file_path) const {
    FILE *fd = fopen(file_path.c_str(), "w");

    fputs("Genomes:\n", fd);
    for (auto it = env_->init_genomes_paths_.begin(); it != env_->init_genomes_paths_.end(); ++it) {
      fputs(it->c_str(), fd);
      fputs("\n", fd);
    }
    fputs("Refining k_sequence:\n", fd);
    for (size_t i = 0; i < env_->k_history_.size(); ++i) {
      if (i != 0) {
        fputs(" -> ", fd);
      }
      fprintf(fd, "%d (#%ld)", env_->k_history_[i], env_->num_genomes_history_[i]);
    }
    fputs("\n", fd);

    fclose(fd);
  }

  void PrepareDirForSave(std::string path) const {
    cap::utils::MakeDirPath(path);
    WriteStateDesc(path + cap_cfg::get().desc_file_name);
  }

  template <class gp_t>
  shared_ptr<gp_t> BuildGPFromStreams(ContigStreams &streams,
                                      unsigned k) const {
    typedef NewExtendedSequenceMapper<Graph, typename gp_t::index_t> Mapper;

    shared_ptr<gp_t> result(new gp_t(k, env_->kDefaultGPWorkdir, 0));

    //fixme use rc_wrapper
    ContigStreams rc_contigs = io::RCWrap(streams);
    rc_contigs.reset();

    debruijn_graph::ConstructGraphUsingExtentionIndex(config::debruijn_config::construction(),
                                                      rc_contigs, result->g, result->index);

    env_->coloring_ = std::make_shared<ColorHandler<Graph> >(result->g, streams.size());
    ColoredGraphConstructor<Graph, Mapper> colored_graph_constructor(result->g,
        *(env_->coloring_), *MapperInstance<gp_t>(*result));
    colored_graph_constructor.ConstructGraph(rc_contigs);

    INFO("Filling positions");
    FillPositions(*result, rc_contigs, env_->coordinates_handler_);
    INFO("Filling positions done.");

    return result;
  }

  //template <class gp_t>
  //shared_ptr<gp_t> BuildGPFromSaves(const size_t K, const std::string &/* path */) const;

  shared_ptr<RtSeqGP> BuildGPFromSaves(const size_t K, const std::string &path) const {
    typedef RtSeqGP gp_t;

    shared_ptr<gp_t> result(new gp_t(unsigned(K), env_->kDefaultGPWorkdir, 0));

    debruijn_graph::graphio::ScanGraphPack(path, *result);

    ContigStreams streams;
    for (size_t i = 0; i < env_->genomes_.size(); ++i) {
      streams.push_back(make_shared<io::SequenceReadStream<Contig>>(
                    env_->genomes_[i], env_->genomes_names_[i]));
    }
    ContigStreams rc_contigs = io::RCWrap(streams);
    rc_contigs.reset();

    INFO("Filling positions");
    FillPositions(*result, rc_contigs, env_->coordinates_handler_);
    INFO("Filling positions done.");

    return result;
  }

  template <class gp_t>
  void SaveCurrentStreams(const gp_t &/* gp */, const std::string &dir) const {
        for (size_t i = 0; i < env_->genomes_.size(); ++i) {
      std::string output_filename = dir + path::filename(env_->init_genomes_paths_[i]);
            if (!output_filename.empty()) {
                Contig contig;
                io::osequencestream out_stream(output_filename);
                DEBUG("Saving to " << output_filename);

        io::SequenceReadStream<io::SingleRead> stream(env_->genomes_[i], env_->genomes_names_[i]);
                while (!stream.eof()) {
                    stream >> contig;
                    out_stream << contig;
                }
            }
        }
  }

  void UpdateStreams() {
    for (unsigned i = 0; i < env_->genomes_.size(); ++i) {
      env_->genomes_[i] = env_->coordinates_handler_.ReconstructGenome(2 * i);
      //VERIFY(env_->genomes_[i]->IsValid());
    }
  }

  template <class gp_t>
  void RefineTemplated(gp_t &gp) {
    INFO("Store threads");
    //env_->coordinates_handler_.StoreGenomeThreads();
    INFO("Store threads ended");
    double delta = 5.;

    //outdated!!!
    //debruijn_config::simplification::bulge_remover br_config;
    //br_config.max_bulge_length_coefficient = 3;
    //br_config.max_coverage = 1000.;
    //br_config.max_relative_coverage = 1.2;
    //br_config.max_delta = delta;
    //br_config.max_relative_delta = 0.1;

    INFO("Removing bulges");

    BulgeRemoverCallbackToCoordinatesHandlerAdapter<Graph> adapter(
        env_->coordinates_handler_);
    boost::function<void(EdgeId, const std::vector<EdgeId> &)> projecting_callback =
        boost::bind(&BulgeRemoverCallbackToCoordinatesHandlerAdapter<Graph>::Project,
                    &adapter, _1, _2);

    //omp_set_num_threads(1);
    debruijn::simplification::RemoveBulges(gp.g, br_config, projecting_callback);
    //omp_set_num_threads(4);

    INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");

    /*
    debruijn_config::simplification::complex_bulge_remover cbr_config;
    cbr_config.enabled = true;
    cbr_config.pics_enabled = false;
    cbr_config.folder = "";
    cbr_config.max_relative_length = 3;
    cbr_config.max_length_difference = 1000;

    INFO("Removing complex bulges");
    RemoveComplexBulges(gp.g, cbr_config);

    INFO("Clipping tips with projection");

    TipsProjector<gp_t> tip_projector(gp);
    boost::function<void(EdgeId)> projecting_callback = boost::bind(
        &TipsProjector<gp_t>::ProjectTip, &tip_projector, _1);
    debruijn_config::simplification::tip_clipper tc_config;

    tc_config.condition = "{ tc_lb 2. }";

    ClipTipsWithProjection(gp, tc_config, true);
    */

    INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");


    env_->k_history_.push_back(env_->GetGraphK());
    env_->num_genomes_history_.push_back(env_->init_genomes_paths_.size());
    env_->coordinates_handler_.DumpRanges();

    UpdateStreams();
  }

  template <class gp_t>
  void FindIndelsTemplated(gp_t& gp, std::ofstream &out_stream,
      const bool mask_indels) {
    GenomeContiguousPathsGraphTraversalConstraints<Graph> traversal_constraints(
        env_->coordinates_handler_);
    SimpleIndelFinder<gp_t> indel_finder(gp, *env_->coloring_,
        env_->coordinates_handler_, traversal_constraints, out_stream,
        mask_indels);
    indel_finder.FindIndelEvents();

    if (mask_indels) {
      env_->k_history_.push_back(env_->GetGraphK());
      env_->num_genomes_history_.push_back(env_->init_genomes_paths_.size());
      env_->coordinates_handler_.DumpRanges();
      UpdateStreams();

    }

  }

  template <class gp_t>
  void FindInversionsTemplated(gp_t& gp, const std::string &base_pic_file_name,
      const bool mask_inversions) const {
    SimpleInversionFinder<gp_t> finder(gp, *env_->coloring_, env_->coordinates_handler_,
        base_pic_file_name, mask_inversions);
    finder.FindInversionEvents();
  }

  template<class gp_t>
  void RefillPositions(const gp_t &gp) {
    ContigStreams streams;
    for (size_t i = 0; i < env_->genomes_.size(); ++i) {
      streams.push_back(make_shared<io::SequenceReadStream<Contig>>(
                    env_->genomes_[i], env_->genomes_names_[i]));
    }
    ContigStreams rc_contigs = io::RCWrap(streams);
    rc_contigs.reset();

    INFO("Filling positions");
    FillPositions(gp, rc_contigs, env_->coordinates_handler_);
    INFO("Filling positions done.");
  }

 public:
  CapEnvironmentManager(CapEnvironment *env)
      : env_(env) {
        //last_streams_used_() {
  }

  virtual ~CapEnvironmentManager() {
    /*
    for (auto it = last_streams_used_.begin(); it != last_streams_used_.end(); ++it) {
      delete *it;
    }
    */
  }

  void ClearEnvironment() const {
    env_->ClearGP();
  }

  std::string GetDirForCurrentState() const {
    std::string env_dir = env_->dir();
    std::stringstream merged_ks_stream;
    for (size_t i = 0; i < env_->k_history_.size(); ++i) {
      if (i != 0) {
        merged_ks_stream << " ";
      }
      merged_ks_stream << env_->k_history_[i] << " " << env_->num_genomes_history_[i];
    }
    std::string cache_dir = "cache_" +
        cap::utils::GenMD5FromFiles(env_->init_genomes_paths_, merged_ks_stream.str());

    return env_dir + "/" + cache_dir + "/";
  }

  void ConstructGraphFromStreams(ContigStreams &streams, unsigned k) {
    ClearEnvironment();
    env_->CheckConsistency();
    //last_streams_used_ = streams;

    VERIFY(env_->gp_rtseq_ == NULL && env_->gp_lseq_ == NULL);
    if (env_->UseLSeqForThisK(k)) {
        VERIFY(false);
//      env_->SetGraphPack(BuildGPFromStreams<LSeqGP>(
//          streams, k));
    } else {
      env_->SetGraphPack(BuildGPFromStreams<RtSeqGP>(
          streams, k));
    }
  }

  void ConstructGraph(unsigned k) {
    ContigStreams streams;
    for (size_t i = 0; i < env_->genomes_.size(); ++i) {
      streams.push_back(make_shared<io::SequenceReadStream<Contig>>(
                    env_->genomes_[i], env_->genomes_names_[i]));
    }

    ConstructGraphFromStreams(streams, k);
  }

  void SaveGraph(std::string folder) const {
    cap::utils::MakeDirPath(folder);
    VERIFY(cap::utils::DirExist(folder));

    std::string filename = folder + "graph";

    // Saving graph
    /*
    debruijn_graph::graphio::PrinterTraits<Graph>::Printer printer(*env_->graph_);
    printer.SaveGraph(filename);
    printer.SaveEdgeSequences(filename);
    printer.SavePositions(filename, *env_->edge_pos_);
  */
    if (env_->LSeqIsUsed()) {
      //PrintGraphPack(filename, env_->gp_lseq_);
    } else {
      debruijn_graph::graphio::PrintGraphPack(filename, *env_->gp_rtseq_);
    }

    // Saving coloring of graph
    cap::SaveColoring(*env_->graph_, *env_->coloring_, filename);
  }

  void DrawPics(std::string folder) const {
    // Saving pics
    cap::utils::MakeDirPath(folder);
    VERIFY(cap::utils::DirExist(folder));

    std::vector<std::string> genomes_names;
    for (const auto &gname : env_->genomes_names_) {
      genomes_names.push_back(gname);
      genomes_names.push_back(gname + "_RC");
    }
    cap::PrintColoredGraphWithColorFilter(*env_->graph_, *env_->coloring_,
        env_->coordinates_handler_, genomes_names, folder);
  }

  void SetGenomes(const std::vector<std::string> &genomes_paths,
                  const std::vector<std::string> &/* genomes_names */) const {
    VERIFY(env_->init_genomes_paths_.size() == 0);

    env_->init_genomes_paths_ = genomes_paths;
  }

  /*
   * Returns true if added successfully
   */
  bool AddGenomeFromFile(const std::string &filename,
                         const std::string &name,
                         bool crop_repeats = false) const {
    if (!CheckFileExists(filename)) {
      return false;
    }

    if (crop_repeats) {
        JunkCroppingWrapper reader(make_shared<io::FileReadStream>(filename));
        io::SingleRead genome;
        reader >> genome;

        if (!genome.IsValid()) {
            return false;
        }

        env_->init_genomes_paths_.push_back(filename);
        env_->genomes_.push_back(genome.sequence());
        env_->genomes_names_.push_back(name);
        env_->coordinates_handler_.StoreGenomeThreadManual(uint(env_->genomes_.size() - 1),
                                                            reader.coordinates_ladder());
    } else {
        io::FileReadStream reader(filename);
        io::SingleRead genome;
        reader >> genome;

        if (!genome.IsValid()) {
          return false;
        }

        env_->init_genomes_paths_.push_back(filename);
        env_->genomes_.push_back(genome.sequence());
        env_->genomes_names_.push_back(name);
    }


    return true;
  }

  void SaveGenomesToDisk(const bool force) const {
    const std::string &dir = GetDirForCurrentState();

    if (!force && cap::utils::DirExist(dir)) {
      return;
    }
    PrepareDirForSave(dir);

    if (env_->LSeqIsUsed()) {
      SaveCurrentStreams(*env_->gp_lseq_, dir);
    } else {
      SaveCurrentStreams(*env_->gp_rtseq_, dir);
    }
  }

  void Refine() {
    env_->CheckConsistency();

    if (env_->LSeqIsUsed()) {
      RefineTemplated(*env_->gp_lseq_);
    } else {
      RefineTemplated(*env_->gp_rtseq_);
    }

    env_->CheckConsistency();
  }

  int FindIndels(const bool mask_indels, const std::string &output_file,
      const std::string &output_mode) {
    std::ios_base::openmode mode;
    if (output_mode == "w") {
      mode = std::ios_base::out;
    } else {
      mode = std::ios_base::app;
    }
    std::ofstream stream(output_file, mode);

    if (stream.fail()) {
      return 1;
    }

    if (env_->LSeqIsUsed()) {
      FindIndelsTemplated(*env_->gp_lseq_, stream, mask_indels);
    } else {
      FindIndelsTemplated(*env_->gp_rtseq_, stream, mask_indels);
    }

    return 0;
  }

  int FindInversions(const bool mask_inversions, const std::string &/* output_file */, const std::string &/* output_mode */) const {
    const std::string &dir = GetDirForCurrentState();
    const std::string &base_pic_dir = dir + "/inversions";
    const std::string &base_pic_file_name = base_pic_dir + "/inv";

    cap::utils::MakeDirPath(base_pic_dir);
    /*
    std::ios_base::openmode mode;
    if (output_mode == "w") {
      mode = std::ios_base::out;
    } else {
      mode = std::ios_base::app;
    }
    std::ofstream stream(output_file, mode);

    if (stream.fail()) {
      return 1;
    }
    */

    if (env_->LSeqIsUsed()) {
      FindInversionsTemplated(*env_->gp_lseq_, base_pic_file_name, mask_inversions);
    } else {
      FindInversionsTemplated(*env_->gp_rtseq_, base_pic_file_name, mask_inversions);
    }

    return 0;
  }

  void LoadGraphFromSaves(const unsigned K, const std::string &path) {
    ClearEnvironment();
    env_->CheckConsistency();
    //last_streams_used_ = streams;

    VERIFY(env_->gp_rtseq_ == NULL && env_->gp_lseq_ == NULL);
    if (env_->UseLSeqForThisK(K)) {
      //env_->SetGraphPack(BuildGPFromSaves<LSeqGP>(K, path));
    } else {
      env_->SetGraphPack(BuildGPFromSaves(K, path));
    }

    env_->coloring_ = std::make_shared<ColorHandler<Graph> >(env_->graph(), env_->genome_cnt());
    INFO("Loading coloring from " << path);
    cap::LoadColoring(*env_->graph_, *env_->element_finder_, *env_->coloring_, path);

    env_->CheckConsistency();
  }
};

}
