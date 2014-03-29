// FIXME: Refactor and turn into stage

//todo rewrite with extended sequence mapper!
template<class Graph, class Index>
class EtalonPairedInfoCounter {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	const Index& index_;
	const KmerMapper<Graph>& kmer_mapper_;
	size_t k_;

	size_t insert_size_;
	size_t read_length_;
	int gap_;
	size_t delta_;

  void AddEtalonInfo(PairedInfoIndexT<Graph>& index, EdgeId e1, EdgeId e2, double d) {
    index.AddPairInfo(e1, e2, d, 1000., 0.);
	}

  void ProcessSequence(const Sequence& sequence, PairedInfoIndexT<Graph>& index)
  {
		int mod_gap = (gap_ + (int) k_ > (int) delta_ ) ? gap_ - (int) delta_ : 0 - (int) k_;
		runtime_k::RtSeq left(k_ +1, sequence);
		left >>= 0;
		for (size_t left_idx = 0;
             left_idx + 2 * (k_ + 1) + mod_gap <= sequence.size();
             ++left_idx) {
			left <<= sequence[left_idx + k_];
			runtime_k::RtSeq left_upd = kmer_mapper_.Substitute(left);
			if (!index_.contains(left_upd)) {
				continue;
			}
			pair<EdgeId, size_t> left_pos = index_.get(left_upd);

			size_t right_idx = left_idx + k_ + 1 + mod_gap;
			runtime_k::RtSeq right(k_ + 1, sequence, right_idx);
			right >>= 0;
			for (;
			     right_idx + k_ + 1 <= left_idx + insert_size_ + delta_ && right_idx + k_ + 1 <= sequence.size();
			     ++right_idx) {
				right <<= sequence[right_idx + k_];
				runtime_k::RtSeq right_upd = kmer_mapper_.Substitute(right);
				if (!index_.contains(right_upd)) {
					continue;
				}
				pair<EdgeId, size_t> right_pos = index_.get(right_upd);

				AddEtalonInfo(index, left_pos.first, right_pos.first,
				              0. + (double) right_idx - (double) left_idx +
				              (double) left_pos.second - (double) right_pos.second);
			}
		}
	}

public:
    EtalonPairedInfoCounter(const Graph& g, const Index& index,
                            const KmerMapper<Graph>& kmer_mapper,
                            size_t insert_size, size_t read_length,
                            size_t delta, size_t k)
            : g_(g),
              index_(index),
              kmer_mapper_(kmer_mapper),
              k_(k),
              insert_size_(insert_size),
              read_length_(read_length),
              gap_((int) (insert_size_ - 2 * read_length_)),
              delta_(delta) {
//		VERIFY(insert_size_ >= 2 * read_length_);
    }

    void FillEtalonPairedInfo(const Sequence& genome,
                              omnigraph::de::PairedInfoIndexT<Graph>& paired_info) {
        ProcessSequence(genome, paired_info);
        ProcessSequence(!genome, paired_info);
    }
};

template<class Graph>
void GetAllDistances(const PairedInfoIndexT<Graph>& paired_index,
                     PairedInfoIndexT<Graph>& result,
                     const GraphDistanceFinder<Graph>& dist_finder) {
    for (auto iter = paired_index.begin(); iter != paired_index.end(); ++iter) {
        EdgeId e1 = iter.first();
        EdgeId e2 = iter.second();
        vector<size_t> forward = dist_finder.GetGraphDistancesLengths(e1, e2);
        for (size_t i = 0; i < forward.size(); ++i)
            result.AddPairInfo(e1, e2, (double) forward[i], -10.0, 0.0, false);
    }
}

template<class Graph>
void GetAllDistances(const Graph& g,
                     const PairedInfoIndexT<Graph>& paired_index,
                     const PairedInfoIndexT<Graph>& clustered_index,
                     const GraphDistanceFinder<Graph>& dist_finder,
                     PairedInfoIndexT<Graph>& result)
{
    typedef typename Graph::EdgeId EdgeId;
    typedef vector<EdgeId> Path;
    for (auto iter = paired_index.begin(); iter != paired_index.end(); ++iter) {
        EdgeId first = iter.first();
        EdgeId second = iter.second();
        const vector<Path>& raw_paths = dist_finder.GetGraphDistances(first, second);
        // adding first edge to every path
        vector<Path> paths;
        for (size_t i = 0; i < raw_paths.size(); ++i) {
            Path path;
            path.push_back(first);
            for (size_t j = 0; j < raw_paths[i].size(); ++j)
                path.push_back(raw_paths[i][j]);
            path.push_back(second);

            paths.push_back(path);
        }
        vector<size_t> path_lengths;
        vector<double> path_weights;
        for (size_t i = 0; i < paths.size(); ++i) {
            size_t len_total = 0 ;
            double weight_total = 0.;
            for (size_t j = 0; j < paths[i].size(); ++j) {
                len_total += g.length(paths[i][j]);
                size_t cur_length = 0;
                for (size_t l = j + 1; l < paths[i].size(); ++l) {
                    cur_length += g.length(paths[i][l - 1]);
                    const de::Histogram& infos = clustered_index.GetEdgePairInfo(paths[i][j], paths[i][l]);
                    for (auto iterator = infos.begin(); iterator != infos.end(); ++iterator) {
                        const Point& info = *iterator;
                        if (info.d == cur_length) {
                            weight_total += info.weight;
                            break;
                        }
                    }
                }
            }
            path_lengths.push_back(len_total - g.length(second));
            path_weights.push_back(weight_total);
        }

        for (size_t i = 0; i < paths.size(); ++i) {
            cout << first.int_id() << "(" << g.length(first) << ") "
                 << second.int_id() << "(" << g.length(second) << ") : "
                 << (i + 1) << "-th path (" << path_lengths[i] << ", " << path_weights[i] << ")   :::   ";
            for (size_t j = 0; j < paths[i].size(); ++j) {
                cout << paths[i][j].int_id() << "(" << g.length(paths[i][j]) << ") ";
            }
            cout << endl;
        }
    }
}

template<class Graph, class Index>
void FillEtalonPairedIndex(PairedInfoIndexT<Graph>& etalon_paired_index,
                           const Graph &g, const Index& index,
                           const KmerMapper<Graph>& kmer_mapper, size_t is, size_t rs,
                           size_t delta, const Sequence& genome, size_t k)
{
    VERIFY_MSG(genome.size() > 0,
               "The genome seems not to be loaded, program will exit");
    INFO((string) (FormattedString("Counting etalon paired info for genome of length=%i, k=%i, is=%i, rs=%i, delta=%i")
                   << genome.size() << k << is << rs << delta));

    EtalonPairedInfoCounter<Graph, Index> etalon_paired_info_counter(g, index, kmer_mapper, is, rs, delta, k);
    etalon_paired_info_counter.FillEtalonPairedInfo(genome, etalon_paired_index);

    DEBUG("Etalon paired info counted");
}

template<class Graph, class Index>
void FillEtalonPairedIndex(PairedInfoIndexT<Graph>& etalon_paired_index,
                           const Graph &g, const Index& index,
                           const KmerMapper<Graph>& kmer_mapper, const Sequence& genome,
                           const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                           size_t k) {

    FillEtalonPairedIndex(etalon_paired_index, g, index, kmer_mapper,
                          size_t(lib.data().mean_insert_size), lib.data().read_length, size_t(lib.data().insert_size_deviation),
                          genome, k);

    //////////////////DEBUG
    //	SimpleSequenceMapper<k + 1, Graph> simple_mapper(g, index);
    //	Path<EdgeId> path = simple_mapper.MapSequence(genome);
    //	SequenceBuilder sequence_builder;
    //	sequence_builder.append(Seq<k>(g.EdgeNucls(path[0])));
    //	for (auto it = path.begin(); it != path.end(); ++it) {
    //		sequence_builder.append(g.EdgeNucls(*it).Subseq(k));
    //	}
    //	Sequence new_genome = sequence_builder.BuildSequence();
    //	NewEtalonPairedInfoCounter<k, Graph> new_etalon_paired_info_counter(g, index,
    //			insert_size, read_length, insert_size * 0.1);
    //	PairedInfoIndexT<Graph> new_paired_info_index(g);
    //	new_etalon_paired_info_counter.FillEtalonPairedInfo(new_genome, new_paired_info_index);
    //	CheckInfoEquality(etalon_paired_index, new_paired_info_index);
    //////////////////DEBUG
    //	INFO("Etalon paired info counted");
}

template<class Graph>
void CountPairedInfoStats(const Graph& g,
                          const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                          const PairedInfoIndexT<Graph>& paired_index,
                          const PairedInfoIndexT<Graph>& etalon_index,
                          const string& output_folder) {
    PairedInfoIndexT<Graph> filtered_index = paired_index;
//    PairInfoWeightFilter<Graph>(g, 40).Filter(filtered_index);
    PairInfoFilter<Graph>(PairInfoWeightChecker<Graph>(g, 40)).Filter(filtered_index);
    INFO("Counting paired info stats");
    EdgePairStat<Graph>(g, paired_index, output_folder).Count();

    //todo remove filtration if launch on etalon info is ok
    UniquePathStat<Graph>(g, filtered_index,
                          (size_t)math::round(lib.data().mean_insert_size),
                          lib.data().read_length,
                          0.1 * lib.data().mean_insert_size).Count();
    UniqueDistanceStat<Graph>(etalon_index).Count();
    INFO("Paired info stats counted");
}

// leave only those pairs, which edges have no path in the graph between them
template<class Graph>
void FilterIndexWithExistingPaths(PairedIndexT& scaf_clustered_index,
                                  const PairedIndexT& index,
                                  const conj_graph_pack &gp,
                                  const GraphDistanceFinder<Graph>& dist_finder) {
    for (auto it = index.begin(); it != index.end(); ++it) {
        const de::Histogram& histogram = *it;
        EdgeId e1 = it.first();
        EdgeId e2 = it.second();
        if (gp.g.OutgoingEdgeCount(gp.g.EdgeEnd(e1)) == 0 && gp.g.IncomingEdgeCount(gp.g.EdgeEnd(e1)) == 1 &&
            gp.g.IncomingEdgeCount(gp.g.EdgeStart(e2)) == 0 && gp.g.OutgoingEdgeCount(gp.g.EdgeStart(e2)) == 1)     {
            vector<size_t> dists = dist_finder.GetGraphDistancesLengths(e1, e2);
            if (dists.size() == 0)
                for (auto point_iter = histogram.begin(); point_iter != histogram.end(); ++point_iter)
                    if (math::gr(point_iter->d, 0.)) {
                        scaf_clustered_index.AddPairInfo(it.first(), it.second(),
                                                         point_iter->d, point_iter->weight, 20.);
                    }
        }
    }
}

inline
void tSeparatedStats(conj_graph_pack& gp, const Sequence& contig,
                     PairedInfoIndex<conj_graph_pack::graph_t> &ind,
                     const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                     size_t /*k*/) {
    typedef omnigraph::de::PairInfo<EdgeId> PairInfo;

    MappingPath<Graph::EdgeId> m_path1 = FindGenomeMappingPath(contig, gp.g,
                                                               gp.index, gp.kmer_mapper);

    map<Graph::EdgeId, vector<pair<int, int>>> inGenomeWay;
    int CurI = 0;
    int gaps = 0;
    for (size_t i = 0; i < m_path1.size(); i++) {
        bool new_edge_added = false;
        EdgeId ei = m_path1[i].first;
        MappingRange mr = m_path1[i].second;
        int start = (int)(mr.initial_range.start_pos - mr.mapped_range.start_pos);
        if (inGenomeWay.find(ei) == inGenomeWay.end()) {
            vector<pair<int, int>> tmp;
            tmp.push_back(make_pair(CurI, start));
            inGenomeWay[ei] = tmp;
            CurI++;
            new_edge_added = true;
            DEBUG("Edge " << gp.g.str(ei) << " num " << CurI << " pos " << start);
        } else {
            if (m_path1[i - 1].first == ei) {
                if (abs(start - inGenomeWay[ei][(inGenomeWay[ei].size() - 1)].second) > 50) {
                    inGenomeWay[ei].push_back(make_pair(CurI, start));
                    CurI++;
                    new_edge_added = true;
                    DEBUG("Edge " << gp.g().str(ei) << " num " << CurI << " pos " << start);
                }
            } else {
                inGenomeWay[ei].push_back(make_pair(CurI, start));
                CurI++;
                new_edge_added = true;
                DEBUG("Edge " << gp.g.str(ei) << " num " << CurI << " pos " << start);
            }
        }
        if (new_edge_added && (i > 0)) {
            if (gp.g.EdgeStart(ei) != gp.g.EdgeEnd(m_path1[i - 1].first)) {
                gaps++;
            }
        }
    }
    INFO("Totaly " << CurI << " edges in genome path, with " << gaps << "not adjacent conequences");

    vector<int> stats(10);
    vector<int> stats_d(10);
    int PosInfo = 0;
    int AllignedPI = 0;
    int ExactDPI = 0;
    int OurD = int(lib.data().mean_insert_size) - int(lib.data().read_length);
    for (auto p_iter = ind.begin(), p_end_iter = ind.end();
         p_iter != p_end_iter; ++p_iter) {
        vector<PairInfo> pi = *p_iter;
        for (size_t j = 0; j < pi.size(); j++) {
            EdgeId left_edge = pi[j].first;
            EdgeId right_edge = pi[j].second;
            double d = pi[j].d();
            if (d < 0.001)
                continue;
            int best_d = 100;
            int best_t = 0;
            PosInfo++;
            DEBUG(
                "PairInfo " << gp.g().str(left_edge) << " -- " << gp.g().str(right_edge) << " d " << d);
            bool ExactOnD = false;
            for (size_t left_i = 0; left_i < inGenomeWay[left_edge].size();
                 left_i++)
                for (size_t right_i = 0;
                     right_i < inGenomeWay[right_edge].size(); right_i++) {
                    if (best_d
                        > abs(
                            inGenomeWay[right_edge][right_i].second
                            - inGenomeWay[left_edge][left_i].second
                            - d)) {
                        best_d = (int)math::round(abs(
                            inGenomeWay[right_edge][right_i].second
                            - inGenomeWay[left_edge][left_i].second
                            - d));
                        best_t = inGenomeWay[right_edge][right_i].first
                                 - inGenomeWay[left_edge][left_i].first;
                        DEBUG("best d " << best_d);
                        if ((inGenomeWay[right_edge][right_i].second
                             - inGenomeWay[left_edge][left_i].second
                             - (int) gp.g.length(left_edge) <= OurD)
                            && (inGenomeWay[right_edge][right_i].second
                                - inGenomeWay[left_edge][left_i].second
                                + (int) gp.g.length(right_edge) >= OurD))
                            ExactOnD = true;
                        else
                            ExactOnD = false;
                    }
                }
            if (best_t > 5)
                best_t = 5;
            if (best_d < 100) {
                AllignedPI++;
                stats[best_t]++;
                if (ExactOnD) {
                    stats_d[best_t]++;
                    ExactDPI++;
                }
            }

        }
    }INFO(
        "Total positive pair info " << PosInfo << " alligned to genome " << AllignedPI << " with exact distance " << ExactDPI);
    INFO(
        "t-separated stats Alligneg: 1 - " << stats[1] << " 2 - " << stats[2] << " 3 - " << stats[3] << " 4 - " << stats[4] << " >4 - " << stats[5]);
    INFO(
        "t-separated stats Exact: 1 - " << stats_d[1] << " 2 - " << stats_d[2] << " 3 - " << stats_d[3] << " 4 - " << stats_d[4] << " >4 - " << stats[5]);
}

template<class Graph>
void CountAndSaveAllPaths(const Graph& g, const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                          const PairedInfoIndexT<Graph>& paired_index, const PairedInfoIndexT<Graph>& /*clustered_index*/) {
    PairedIndexT all_paths(g);
    GetAllDistances<Graph>(paired_index,
                           all_paths,
                           GraphDistanceFinder<Graph>(g,
                                                      size_t(lib.data().mean_insert_size),
                                                      lib.data().read_length,
                                                      size_t(lib.data().insert_size_deviation)));

    std::string dir_name = cfg::get().output_dir + "estimation_qual/";
    make_dir(dir_name);

    typename PrinterTraits<Graph>::Printer printer(g);
    printer.savePaired(dir_name + "paths", all_paths);

    //PairedIndexT& all_paths_2(g);
    //GetAllDistances<Graph>(g,
    //paired_index, clustered_index,
    //all_paths_2,
    //GraphDistanceFinder<Graph>(g, *cfg::get().ds.IS, *cfg::get().ds.RL,
    //size_t(*cfg::get().ds.is_var)));
    //printer.savePaired(dir_name + "paths_all", all_paths_2);
}

void FillAndCorrectEtalonPairedInfo(PairedIndexT&  corrected_etalon_index,
                                    const conj_graph_pack& gp,
                                    const PairedIndexT&  paired_index, size_t insert_size,
                                    size_t read_length, size_t delta,
                                    bool save_etalon_info_history = false) {
    INFO("Filling etalon paired index");
    PairedIndexT etalon_index(gp.g);
    bool successful_load = false;
    if (cfg::get().entry_point >= ws_distance_estimation) {
        string p = path::append_path(cfg::get().load_from, "../etalon");
        if (!path::is_regular_file(p + ".prd")) {
            DEBUG("file " << p + ".prd" << " does not exist");
        }
        else {
            INFO("Loading etalon pair info from the previous run...");
            Graph& graph = const_cast<Graph&>(gp.g);
            ScannerTraits<Graph>::Scanner scanner(graph);
            scanner.loadPaired(p, etalon_index);
            path::files_t files;
            files.push_back(p);
            path::copy_files_by_prefix(files, cfg::get().output_dir);
            successful_load = true;
        }
    }
    if (!successful_load)
        FillEtalonPairedIndex(etalon_index, gp.g,
                              gp.index, gp.kmer_mapper, insert_size, read_length, delta,
                              gp.genome, gp.k_value);
    INFO("Etalon paired index filled");

    INFO("Correction of etalon paired info has been started");

    INFO("Filtering etalon info");
    //leave only info between edges both present in paired_index
    PairedIndexT filtered_etalon_index(gp.g);
    for (auto iter = etalon_index.begin(); iter != etalon_index.end(); ++iter) {
        const de::Histogram& histogram = *iter;
        EdgeId first_edge = iter.first();
        EdgeId second_edge = iter.second();
        if (paired_index.GetEdgePairInfo(first_edge, second_edge).size() > 0) {
            for (auto point = histogram.begin(); point != histogram.end(); ++point)
                filtered_etalon_index.AddPairInfo(first_edge, second_edge, *point);
        }
        else
            DEBUG("Filtering out pair_info " << gp.g.int_id(first_edge) << " "
                  << gp.g.int_id(second_edge));
    }

    INFO("Pushing etalon info through estimator");
    GraphDistanceFinder<Graph> dist_finder(gp.g, insert_size, read_length, delta);
    DistanceEstimator<Graph> estimator(gp.g, filtered_etalon_index, dist_finder, 0., 4.);
    estimator.Estimate(corrected_etalon_index);
    if (save_etalon_info_history) {
        INFO("Saving etalon paired info indices on different stages");
        ConjugateDataPrinter<Graph> data_printer(gp.g);
        data_printer.savePaired(cfg::get().output_dir + "etalon", etalon_index);
        data_printer.savePaired(cfg::get().output_dir + "etalon_filtered_by_index",
                                filtered_etalon_index);
        data_printer.savePaired(cfg::get().output_dir + "etalon_corrected_by_graph",
                                corrected_etalon_index);
        INFO("Everything is saved");

        if (cfg::get().paired_info_scaffolder) {
            GraphDistanceFinder<Graph> dist_finder(gp.g, insert_size, read_length, delta);
            INFO("Saving paired information statistics for a scaffolding");
            PairedIndexT scaf_etalon_index(gp.g);
            FilterIndexWithExistingPaths(scaf_etalon_index, etalon_index, gp, dist_finder);
            data_printer.savePaired(
                cfg::get().output_dir + "scaf_etalon",
                scaf_etalon_index);
            PairedIndexT scaf_filtered_etalon_index(gp.g);
            FilterIndexWithExistingPaths(scaf_filtered_etalon_index, filtered_etalon_index, gp, dist_finder);
            data_printer.savePaired(
                cfg::get().output_dir + "scaf_etalon_filtered",
                scaf_filtered_etalon_index);
        }

        INFO("Everything saved");
    }
    INFO("Correction finished");
}

void CountClusteredPairedInfoStats(const conj_graph_pack &gp,
                                   const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                                   const PairedInfoIndexT<Graph> &paired_index,
                                   const PairedInfoIndexT<Graph> &clustered_index) {
    PairedIndexT etalon_index(gp.g);

    FillAndCorrectEtalonPairedInfo(etalon_index, gp, paired_index,
                                 (size_t)math::round(lib.data().mean_insert_size),
                                 lib.data().read_length,
                                 (size_t)math::round(lib.data().insert_size_deviation), true);

	CountAndSaveAllPaths(gp.g, lib, paired_index, clustered_index);

	INFO("Counting clustered info stats");
	EdgeQuality<Graph, Index> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
  //EstimationQualityStat<Graph> estimation_stat(gp.g, edge_qual,
                                              //paired_index, clustered_index, etalon_index);
  //estimation_stat.Count();
  //estimation_stat.SaveStats(cfg::get().output_dir + "estimation_qual/");

	INFO("Counting overall cluster stat");
	ClusterStat<Graph>(clustered_index).Count();
	INFO("Overall cluster stat");

  if (cfg::get().paired_info_scaffolder) {
		ConjugateDataPrinter<Graph> data_printer(gp.g);
    INFO("Generating the statistics of pair info for scaffolding");
    PairedIndexT scaf_clustered_index(gp.g);
    FilterIndexWithExistingPaths(scaf_clustered_index,
                                 clustered_index, gp,
                                 GraphDistanceFinder<Graph>(gp.g,
                                         (size_t)math::round(lib.data().mean_insert_size),
                                         lib.data().read_length,
                                         (size_t)math::round(lib.data().insert_size_deviation)));
    data_printer.savePaired(cfg::get().output_dir + "scaf_clustered",
                            scaf_clustered_index);
  }
  //  PairedInfoIndexT<Graph> etalon_clustered_index;
	//	DistanceEstimator<Graph> estimator(g, etalon_index, insert_size,
	//			max_read_length, cfg::get().de.delta,
	//			cfg::get().de.linkage_distance, cfg::get().de.max_distance);
	//	estimator.Estimate(etalon_clustered_index);

  //  PairedInfoIndexT<Graph> filtered_clustered_index(g);
	//	PairInfoFilter<Graph> (g, 1000.).Filter(
  //      clustered_index[>etalon_clustered_index<], filtered_clustered_index);
	INFO("Counting mate-pair transformation stat");
	MatePairTransformStat<Graph>(gp.g, //filtered_
	    clustered_index).Count();
	INFO("Mate-pair transformation stat counted");
	INFO("Clustered info stats counted");
}
