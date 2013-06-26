#pragma once

#include "graph_simplification.hpp"

#include "compare_standard.hpp"
#include "graphio.hpp"

#include "comparison_utils.hpp"
#include "diff_masking.hpp"
#include "repeat_masking.hpp"
#include "genome_correction.hpp"
#include "assembly_compare.hpp"
#include "simple_indel_finder.hpp"
#include "simple_inversion_finder.hpp"
#include "test_utils.hpp"

#include "cap_environment.hpp"
#include "io/sequence_reader.hpp"
#include "omni/loop_killer.hpp"

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
                                      unsigned k, bool fill_pos = true) const {
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

    for (auto ptr : rc_contigs) {
      delete ptr;
    }

    INFO("Filling positions");
    if (fill_pos) {
      env_->coordinates_handler_.SetGraph(&(result->g));

      unsigned char genome_id = 0;

      for (auto it = streams.begin(); it != streams.end(); ++it) {
        cap::RCWrapper stream(**it);
        stream.reset();

        Mapper mapper(result->g, result->index, result->kmer_mapper, result->k_value + 1);
        io::SingleRead genome;
        // for forward and reverse directions
        for (int i = 0; i < 2; ++i) {
          stream >> genome;
          MappingPath<EdgeId> mapping_path = mapper.MapRead(genome);
          const std::vector<EdgeId> edge_path =
              mapping_path.simple_path().sequence();
          env_->coordinates_handler_.AddGenomePath(2 * genome_id + i, edge_path);
        }

        genome_id++;
      }
    }
    INFO("Filling positions done.");

    return result;
  }

  template <class gp_t>
  shared_ptr<gp_t> BuildGPFromSaves(const size_t K, const std::string path) const {
    shared_ptr<gp_t> result(new gp_t(K, env_->kDefaultGPWorkdir));

    //ScanGraphPack(path, *result);
    // TODO

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

        io::SequenceReader<io::SingleRead> stream(*env_->genomes_[i], env_->genomes_names_[i]);
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
      env_->genomes_[i] = std::make_shared<Sequence>(
          env_->coordinates_handler_.ReconstructGenome(2 * i));
      //VERIFY(env_->genomes_[i]->IsValid());
		}
  }

  template <class gp_t>
  void RefineTemplated(gp_t &gp) {
    INFO("Store threads");
    //env_->coordinates_handler_.StoreGenomeThreads();
    INFO("Store threads ended");
    size_t delta = 5;

    debruijn_config::simplification::bulge_remover br_config;
    br_config.max_bulge_length_coefficient = 3;
    br_config.max_coverage = 1000.;
    br_config.max_relative_coverage = 1.2;
    br_config.max_delta = delta;
    br_config.max_relative_delta = 0.1;

    INFO("Removing bulges");

    BulgeRemoverCallbackToCoordinatesHandlerAdapter<Graph> adapter(
        env_->coordinates_handler_);
    boost::function<void(EdgeId, const std::vector<EdgeId> &)> projecting_callback =
        boost::bind(&BulgeRemoverCallbackToCoordinatesHandlerAdapter<Graph>::Project,
                    &adapter, _1, _2);

    //omp_set_num_threads(1);
    RemoveBulges(gp.g, br_config, projecting_callback);
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

    //INFO("Killing loops");

    //SimpleLoopKiller<typename gp_t::graph_t> loop_killer(gp.g, /*splitting_edge_len*/ env_->GetGraphK() * 3, /*max_component_size*/ 10);
    //loop_killer.KillAllLoops();

    INFO("Remapped " << gp.kmer_mapper.size() << " k-mers");


    env_->k_history_.push_back(env_->GetGraphK());
    env_->num_genomes_history_.push_back(env_->init_genomes_paths_.size());
    env_->coordinates_handler_.StoreGenomeThreads();

    UpdateStreams(gp);
  }

  template <class gp_t>
  void FindIndelsTemplated(gp_t& gp, std::ofstream &out_stream,
      const bool mask_indels) {
    SimpleIndelFinder<gp_t> indel_finder(gp, *env_->coloring_, 
        env_->coordinates_handler_, out_stream, mask_indels);
    indel_finder.FindIndelEvents();

    if (mask_indels) {
      env_->coordinates_handler_.StoreGenomeThreads();
      UpdateStreams(gp);
    }

    //SimpleInDelCorrector<Graph> corrector(gp.g, *env_->coloring_,
    //    (*MapperInstance(gp)).MapSequence(*env_->genomes_[0]).simple_path().sequence(), /*genome_color*/
    //    kRedColorSet, /*assembly_color*/kBlueColorSet);
    //corrector.Analyze();
  }

  template <class gp_t>
  void FindInversionsTemplated(gp_t& gp, const std::string &base_pic_file_name,
      const bool mask_inversions) const {
    SimpleInversionFinder<gp_t> finder(gp, *env_->coloring_, env_->coordinates_handler_,
        base_pic_file_name, mask_inversions);
    finder.FindInversionEvents();
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

  void ConstructGraphFromStreams(const std::vector<ContigStream *> &streams, unsigned k, bool fill_pos) {
    ClearEnvironment();
    env_->CheckConsistency();
    //last_streams_used_ = streams;

    VERIFY(env_->gp_rtseq_ == NULL && env_->gp_lseq_ == NULL);
    if (env_->UseLSeqForThisK(k)) {
      env_->SetGraphPack(BuildGPFromStreams<LSeqGP>(
          streams, k, fill_pos));
    } else {
      env_->SetGraphPack(BuildGPFromStreams<RtSeqGP>(
          streams, k, fill_pos));
    }
  }

  void ConstructGraph(unsigned k, bool fill_pos) {
    std::vector<ContigStream *> streams;
    for (size_t i = 0; i < env_->genomes_.size(); ++i) {
      streams.push_back(new io::SequenceReader<io::SingleRead>(
                    *env_->genomes_[i], env_->genomes_names_[i]));
    }

    ConstructGraphFromStreams(streams, k, fill_pos);
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
    std::vector<std::string> genomes_names;
    for (const auto &gname : env_->genomes_names_) {
      genomes_names.push_back(gname);
      genomes_names.push_back(gname + "_RC");
    }
    cap::PrintColoredGraphWithColorFilter(*env_->graph_, *env_->coloring_,
        env_->coordinates_handler_, genomes_names, folder + "/pics/graph.dot");
  }

  void SetGenomes(const std::vector<std::string> &genomes_paths,
                  const std::vector<std::string> &genomes_names) const {
    VERIFY(env_->init_genomes_paths_.size() == 0);

    env_->init_genomes_paths_ = genomes_paths;
  }

  /*
   * Returns true if added successfully
   */
  bool AddGenomeFromFile(const std::string &filename,
                         const std::string &name) const {
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
    env_->genomes_names_.push_back(name);

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

  int FindInversions(const bool mask_inversions, const std::string &output_file, const std::string &output_mode) const {
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

  void LoadGraphFromSaves(const size_t K, const std::string &path) {
    ClearEnvironment();
    env_->CheckConsistency();
    //last_streams_used_ = streams;

    VERIFY(env_->gp_rtseq_ == NULL && env_->gp_lseq_ == NULL);
    if (env_->UseLSeqForThisK(K)) {
      env_->SetGraphPack(BuildGPFromSaves<LSeqGP>(K, path));
    } else {
      env_->SetGraphPack(BuildGPFromSaves<RtSeqGP>(K, path));
    }

    cap::SaveColoring(*env_->graph_, *env_->int_ids_, *env_->coloring_, path);

    env_->CheckConsistency();
  }
};

}
