#pragma once

#include "graph_simplification.hpp"

#include "compare_standard.hpp"
#include "graphio.hpp"

#include "comparison_utils.hpp"
#include "diff_masking.hpp"
#include "repeat_masking.hpp"
#include "genome_correction.hpp"
#include "assembly_compare.hpp"
#include "test_utils.hpp"

#include "cap_environment.hpp"
#include "io/sequence_reader.hpp"

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
      fprintf(fd, "%d (#%d)", env_->k_history_[i], env_->num_genomes_history_[i]);
    }
    fputs("\n", fd);

    fclose(fd);
  }

  void PrepareDirForSave(std::string path) const {
    cap::utils::MakeDirPath(path);
    WriteStateDesc(path + cap_cfg::get().desc_file_name);
  }

  template <class gp_t>
  shared_ptr<gp_t> BuildGPFromStreams(std::vector<ContigStream *> streams,
                                                        unsigned k) const {
    typedef NewExtendedSequenceMapper<Graph, typename gp_t::seq_t> Mapper;

    shared_ptr<gp_t> result(new gp_t(k, env_->kDefaultGPWorkdir));

    vector<ContigStream*> rc_contigs;
    for (auto it = streams.begin(); it != streams.end(); ++it) {
      rc_contigs.push_back(new RCWrapper(**it));
      rc_contigs.back()->reset();
    }

    io::ReadStreamVector<ContigStream> rc_read_stream_vector(rc_contigs, false);

    debruijn_graph::ConstructGraph(result->k_value, rc_read_stream_vector, result->g, result->index);

    env_->coloring_ = std::make_shared<ColorHandler<Graph> >(result->g, streams.size());
    ColoredGraphConstructor<Graph, Mapper> colored_graph_constructor(result->g,
        *(env_->coloring_), *MapperInstance<gp_t>(*result));
    colored_graph_constructor.ConstructGraph(rc_contigs);

    return result;
  }

  template <class gp_t>
  void SaveCurrentStreams(const gp_t &gp, const std::string &dir) const {
		for (size_t i = 0; i < env_->genomes_.size(); ++i) {
      std::string output_filename = dir + path::filename(env_->init_genomes_paths_[i]);
			if (!output_filename.empty()) {
				Contig contig;
				io::ofastastream out_stream(output_filename);
				DEBUG("Saving to " << output_filename);

        io::SequenceReader<io::SingleRead> stream(*env_->genomes_[i]);
				while (!stream.eof()) {
					stream >> contig;
					out_stream << contig;
				}
			}
		}
  }

  template <class gp_t>
  void UpdateStreams(const gp_t &gp) {
		for (size_t i = 0; i < env_->genomes_.size(); ++i) {
      io::SequenceReader<io::SingleRead> stream(*env_->genomes_[i]);

      io::ModifyingWrapper<io::SingleRead> refined_stream(
          stream, GraphReadCorrectorInstance(gp.g, *MapperInstance(gp)));
      io::SingleRead refined_genome;

      refined_stream >> refined_genome;

      VERIFY(refined_genome.IsValid());

      env_->genomes_[i] = std::make_shared<Sequence>(refined_genome.sequence());
      // I dont like storing and duplicating of strings
      // last_streams_used_[i] = new io::Reader(env_->genomes_[i]);
		}
  }

  template <class gp_t>
  void RefineTemplated(gp_t &gp) {
    // TODO own config?
    size_t delta = 5;

    debruijn_config::simplification::bulge_remover br_config;
    br_config.max_bulge_length_coefficient = 3;
    br_config.max_coverage = 1000.;
    br_config.max_relative_coverage = 1.2;
    br_config.max_delta = delta;
    br_config.max_relative_delta = 0.1;

    INFO("Removing bulges");
    RemoveBulges(gp.g, br_config);

    INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");

    debruijn_config::simplification::complex_bulge_remover cbr_config;
    cbr_config.enabled = true;
    cbr_config.pics_enabled = false;
    cbr_config.folder = "";
    cbr_config.max_relative_length = 3;
    cbr_config.max_length_difference = 1000;

    INFO("Removing complex bulges");
    RemoveComplexBulges(gp.g, cbr_config);

    TipsProjector<gp_t> tip_projector(gp);
    boost::function<void(EdgeId)> projecting_callback = boost::bind(
        &TipsProjector<gp_t>::ProjectTip, &tip_projector, _1);
    debruijn_config::simplification::tip_clipper tc_config;

    tc_config.condition = "{ tc_lb 2. }";

    INFO("Clipping tips with projection");

    ClipTipsWithProjection(gp, tc_config, true);

    INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");


    env_->k_history_.push_back(env_->GetGraphK());
    env_->num_genomes_history_.push_back(env_->init_genomes_paths_.size());

    UpdateStreams(gp);
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
    // shared_ptr deletes automatically
    env_->gp_rtseq_.reset();
    env_->gp_lseq_.reset();
    env_->graph_ = NULL;
    env_->edge_pos_ = NULL;
    env_->int_ids_ = NULL;
    env_->coloring_.reset();
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

  void ConstructGraphFromStreams(const std::vector<ContigStream *> &streams, unsigned k) {
    ClearEnvironment();
    env_->CheckConsistency();
    //last_streams_used_ = streams;

    VERIFY(env_->gp_rtseq_ == NULL && env_->gp_lseq_ == NULL);
    if (env_->UseLSeqForThisK(k)) {
      env_->gp_lseq_ = BuildGPFromStreams<LSeqGP>(
          streams, k);
      env_->graph_ = &(env_->gp_lseq_->g);
      env_->edge_pos_ = &(env_->gp_lseq_->edge_pos);
      env_->int_ids_ = &(env_->gp_lseq_->int_ids);
    } else {
      env_->gp_rtseq_ = BuildGPFromStreams<RtSeqGP>(
          streams, k);
      env_->graph_ = &(env_->gp_rtseq_->g);
      env_->edge_pos_ = &(env_->gp_rtseq_->edge_pos);
      env_->int_ids_ = &(env_->gp_rtseq_->int_ids);
    }

    env_->CheckConsistency();
  }

  void ConstructGraph(unsigned k) {
    std::vector<ContigStream *> streams;
    for (auto it = env_->genomes_.begin(); it != env_->genomes_.end(); ++it) {
      streams.push_back(new io::SequenceReader<Contig>(**it));
    }

    ConstructGraphFromStreams(streams, k);
  }

  void SaveGraph(std::string folder) const {
    cap::utils::MakeDirPath(folder);
    cap::utils::MakeDirPath(folder + "/saves");
    cap::utils::MakeDirPath(folder + "/pics");

    VERIFY(cap::utils::DirExist(folder));
    VERIFY(cap::utils::DirExist(folder + "/saves"));
    VERIFY(cap::utils::DirExist(folder + "/pics"));

    std::string filename = folder + "/saves/graph";

    // Saving graph
		PrinterTraits<Graph>::Printer printer(*env_->graph_, *env_->int_ids_);
		printer.saveGraph(filename);
		printer.saveEdgeSequences(filename);
		printer.savePositions(filename, *env_->edge_pos_);

    // Saving coloring of graph
    cap::SaveColoring(*env_->graph_, *env_->int_ids_, *env_->coloring_, filename);
    // Saving pics
    cap::PrintColoredGraphWithColorFilter(*env_->graph_, *env_->coloring_,
        *env_->edge_pos_, folder + "/pics/graph.dot");
  }

  void SetGenomes(const std::vector<std::string> &genomes_paths) const {
    VERIFY(env_->init_genomes_paths_.size() == 0);

    env_->init_genomes_paths_ = genomes_paths;
  }

  /*
   * Returns true if added successfully
   */
  bool AddGenomeFromFile(const std::string &filename) const {
    if (!CheckFileExists(filename)) {
      return false;
    }

    io::Reader reader(filename);
    io::SingleRead genome;
    reader >> genome;

    if (!genome.IsValid()) {
      return false;
    }

    env_->init_genomes_paths_.push_back(filename);
    env_->genomes_.push_back(std::make_shared<Sequence>(genome.sequence()));

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


};

}
