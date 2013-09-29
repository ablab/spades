//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "dataset_readers.hpp"
#include "read_converter.hpp"

#include "de/insert_size_refiner.hpp"
#include "de/paired_info.hpp"

#include "utils.hpp"

#include "pair_info_count.hpp"

namespace debruijn_graph {
typedef io::ReadStreamVector<io::SequencePairedReadStream> MultiStreamType;
typedef io::ReadStreamVector<io::PairedReadStream> SingleStreamType;

/**
 * As for now it ignores sophisticated case of repeated consecutive
 * occurrence of edge in path due to gaps in mapping
 *
 * todo talk with Anton about simplification and speed-up of procedure with little quality loss
 */
template<class Graph, class SequenceMapper, class PairedStream>
class LatePairedIndexFiller {
    typedef typename Graph::EdgeId EdgeId;
    typedef runtime_k::RtSeq Kmer;
    typedef boost::function<double(MappingRange, MappingRange)> WeightF;

  public:
    LatePairedIndexFiller(const Graph &graph, const SequenceMapper& mapper, PairedStream& stream, WeightF weight_f) :
            graph_(graph), mapper_(mapper), streams_(stream), weight_f_(weight_f)
    {}

    LatePairedIndexFiller(const Graph &graph, const SequenceMapper& mapper, io::ReadStreamVector< PairedStream >& streams, WeightF weight_f) :
            graph_(graph), mapper_(mapper), streams_(streams), weight_f_(weight_f)
    {}

    bool FillIndex(omnigraph::de::PairedInfoIndexT<Graph>& paired_index) {
        return (streams_.size() == 1 ? FillUsualIndex(paired_index) :  FillParallelIndex(paired_index));
    }

  private:
    template<class PairedRead>
    void ProcessPairedRead(omnigraph::de::PairedInfoIndexT<Graph>& paired_index, const PairedRead& p_r) {
        Sequence read1 = p_r.first().sequence();
        Sequence read2 = p_r.second().sequence();

        MappingPath<EdgeId> path1 = mapper_.MapSequence(read1);
        MappingPath<EdgeId> path2 = mapper_.MapSequence(read2);
        size_t read_distance = p_r.distance();
        for (size_t i = 0; i < path1.size(); ++i) {
            std::pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
            for (size_t j = 0; j < path2.size(); ++j) {
                std::pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
                double weight = weight_f_(mapping_edge_1.second,
                                          mapping_edge_2.second);
                size_t kmer_distance = read_distance
                                       + mapping_edge_2.second.initial_range.end_pos
                                       - mapping_edge_1.second.initial_range.start_pos;
                int edge_distance = (int) kmer_distance
                                    + (int) mapping_edge_1.second.mapped_range.start_pos
                                    - (int) mapping_edge_2.second.mapped_range.end_pos;

                paired_index.AddPairInfo(mapping_edge_1.first,
                                         mapping_edge_2.first,
                                         (double) edge_distance, weight, 0.);
            }
        }
    }

    /**
     * Method reads paired data from stream, maps it to genome and stores it in this PairInfoIndex.
     */
    bool FillUsualIndex(omnigraph::de::PairedInfoIndexT<Graph>& paired_index) {
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it)
            paired_index.AddPairInfo(*it, *it, 0., 0., 0.);

        size_t initial_size = paired_index.size();
        INFO("Processing paired reads (takes a while)");
        PairedStream& stream = streams_.back();
        stream.reset();
        size_t n = 0;
        while (!stream.eof()) {
            typename PairedStream::read_type p_r;
            stream >> p_r;
            ProcessPairedRead(paired_index, p_r);
            VERBOSE_POWER(++n, " paired reads processed");
        }

        return paired_index.size() > initial_size;
    }

    bool FillParallelIndex(omnigraph::de::PairedInfoIndexT<Graph>& paired_index) {
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it)
            paired_index.AddPairInfo(*it, *it, 0., 0., 0.);
        size_t initial_size = paired_index.size();

        INFO("Processing paired reads (takes a while)");
        size_t nthreads = streams_.size();
        std::vector<omnigraph::de::PairedInfoIndexT<Graph>*> buffer_pi(nthreads);

        for (size_t i = 0; i < nthreads; ++i)
            buffer_pi[i] = new omnigraph::de::PairedInfoIndexT<Graph>(graph_);

        size_t counter = 0;
        static const double coeff = 1.3;
#pragma omp parallel num_threads(nthreads)
        {
#pragma omp for reduction(+ : counter)
            for (size_t i = 0; i < nthreads; ++i)
            {
                size_t size = 0;
                size_t limit = 1000000;
                typename PairedStream::read_type r;
                PairedStream& stream = streams_[i];
                stream.reset();
                bool end_of_stream = false;

                //DEBUG("Starting " << omp_get_thread_num());
                while (!end_of_stream) {
                    end_of_stream = stream.eof();
                    while (!end_of_stream && size < limit) {
                        stream >> r;
                        ++counter;
                        ++size;
                        //DEBUG("Processing paired read " << omp_get_thread_num());
                        ProcessPairedRead(*(buffer_pi[i]), r);
                        end_of_stream = stream.eof();
                    }

#pragma omp critical
                    {
                        DEBUG("Merging " << omp_get_thread_num());
                        paired_index.AddAll(*(buffer_pi[i]));
                        DEBUG("Thread number " << omp_get_thread_num()
                              << " is going to increase its limit by " << coeff
                              << " times, current limit is " << limit);
                    }
                    buffer_pi[i]->Clear();
                    limit = (size_t) (coeff * (double) limit);
                }
            }
            DEBUG("Thread number " << omp_get_thread_num() << " finished");
        }
        INFO("Used " << counter << " paired reads");

        for (size_t i = 0; i < nthreads; ++i)
            DEBUG("Size of " << i << "-th map is " << buffer_pi[i]->size());

        for (size_t i = 0; i < nthreads; ++i)
            delete buffer_pi[i];
        INFO("Index built");

        DEBUG("Size of map is " << paired_index.size());

        return paired_index.size() > initial_size;
    }

  private:
    const Graph& graph_;
    const SequenceMapper& mapper_;
    io::ReadStreamVector<PairedStream>& streams_;
    WeightF weight_f_;

    DECL_LOGGER("LatePairedIndexFiller");
};

template<class PairedRead, class Graph, class Mapper>
bool FillPairedIndexWithReadCountMetric(const Graph &g,
                                        const Mapper& mapper,
                                        PairedInfoIndexT<Graph>& paired_info_index,
                                        io::ReadStreamVector<io::IReader<PairedRead> >& streams) {

    INFO("Counting paired info with read count weight");
    LatePairedIndexFiller<Graph, Mapper, io::IReader<PairedRead>>
            pif(g, mapper, streams, PairedReadCountWeight);

    bool res = pif.FillIndex(paired_info_index);
    DEBUG("Paired info with read count weight counted");
    return res;
}

template<class PairedRead, class Graph, class Mapper>
bool FillPairedIndexWithProductMetric(const Graph &g,
                                      const Mapper& mapper,
                                      PairedInfoIndexT<Graph>& paired_info_index,
                                      io::ReadStreamVector<io::IReader<PairedRead> >& streams) {

    INFO("Counting paired info with product weight");

    LatePairedIndexFiller<Graph, Mapper, io::IReader<PairedRead> >
            pif(g, mapper, streams, KmerCountProductWeight);
    bool res = pif.FillIndex(paired_info_index);
    DEBUG("Paired info with product weight counted");
    return res;
}

template<class graph_pack, class PairedRead, class ConfigType>
bool RefineInsertSize(const graph_pack& gp,
                      io::ReadStreamVector<io::IReader<PairedRead> >& streams,
                      ConfigType& config,
                      size_t edge_length_threshold) {
  size_t rl;
  double mean;
  double delta;
  double median;
  double mad;
  std::map<size_t, size_t> percentiles;
  std::map<int, size_t> hist;
  // calling default method
  refine_insert_size(streams, gp, edge_length_threshold, rl, mean, delta, median, mad, percentiles, hist);

  if (hist.size() == 0) {
    config.paired_mode = false;
    WARN("Failed to estimate the insert size of paired reads, because none of the paired reads aligned to long edges.");
    WARN("Paired reads will not be used.");
    return false;
  }

  config.ds.set_IS(mean);
  config.ds.set_is_var(delta);
  config.ds.set_median(median);
  config.ds.set_mad(mad);
  config.ds.set_hist(hist);
  INFO("Mean Insert Size = " << mean);
  INFO("Insert Size stddev= " << delta);
  INFO("Median Insert Size = " << median);
  INFO("Insert Size MAD = " << mad);
  DEBUG("Delta_Mad = " << 1.4826 * mad);

  return true;
}

template<class graph_pack, class PairedRead, class DataSet>
bool RefineInsertSizeForLib(const graph_pack& gp,
                      io::ReadStreamVector<io::IReader<PairedRead> >& streams,
                      DataSet& data,
                      size_t edge_length_threshold) {

  std::map<size_t, size_t> percentiles;
  // calling default method
  data.read_length = 0;
  refine_insert_size(streams, gp, edge_length_threshold,
          data.read_length,
          data.mean_insert_size,
          data.insert_size_deviation,
          data.median_insert_size,
          data.insert_size_mad,
          percentiles,
          data.insert_size_distribution);

  if (data.insert_size_distribution.size() == 0) {
    return false;
  }

  return true;
}

void PairInfoCount::run(conj_graph_pack &gp) {
    if (!cfg::get().developer_mode) {
        gp.paired_indices.Attach();
        gp.paired_indices.Init();
    }

    size_t edge_length_threshold = Nx(gp.g, 50);
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd ||
            cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs) {

            bool insert_size_refined;
            if (cfg::get().use_multithreading) {
                auto streams = paired_binary_readers(cfg::get().ds.reads[i], false, 0);
                insert_size_refined = RefineInsertSizeForLib(gp, *streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
            } else {
                auto_ptr<PairedReadStream> stream = paired_easy_reader(cfg::get().ds.reads[i], false, 0);
                SingleStreamType streams(stream.get());
                streams.release();
                insert_size_refined = RefineInsertSizeForLib(gp, streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
            }

            if (!insert_size_refined) {
                cfg::get_writable().ds.reads[i].data().mean_insert_size = 0.0;

                WARN("Unable to estimate insert size for paired library #" << i);
                if (cfg::get().ds.reads[i].data().read_length <= cfg::get().K) {
                    WARN("Maximum read length (" << cfg::get().ds.reads[i].data().read_length << ") should be greater than K (" << cfg::get().K << ")");
                } else if (cfg::get().ds.reads[i].data().read_length <= cfg::get().K * 11 / 10) {
                    WARN("Maximum read length (" << cfg::get().ds.reads[i].data().read_length << ") is probably too close to K (" << cfg::get().K << ")");
                } else {
                    WARN("None of paired reads aligned properly. Please, check orientation of your read pairs.");
                }

                continue;
            } else {
                INFO("Estimated insert size for paired library #" << i);
                INFO("Insert size = " << cfg::get().ds.reads[i].data().mean_insert_size << ", deviation = " << cfg::get().ds.reads[i].data().insert_size_deviation);
                INFO("Read length = " << cfg::get().ds.reads[i].data().read_length);
            }

            if (cfg::get().use_multithreading) {
                auto paired_streams = paired_binary_readers(cfg::get().ds.reads[i], true, (size_t) cfg::get().ds.reads[i].data().mean_insert_size);
                FillPairedIndexWithReadCountMetric(gp.g, *MapperInstance(gp), gp.paired_indices[i], *paired_streams);
            } else {
                auto_ptr<PairedReadStream> paired_stream = paired_easy_reader(cfg::get().ds.reads[i], true, (size_t) cfg::get().ds.reads[i].data().mean_insert_size);
                SingleStreamType paired_streams(paired_stream.get());
                paired_stream.release();
                FillPairedIndexWithReadCountMetric(gp.g, *MapperInstance(gp), gp.paired_indices[i], paired_streams);
            }
        }
    }
}

}
