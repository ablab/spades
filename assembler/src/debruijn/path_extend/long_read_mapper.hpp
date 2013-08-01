/*
 * long_read_mapper.hpp
 *
 *  Created on: Jun 17, 2013
 *      Author: andrey
 */

#ifndef LONG_READ_MAPPER_HPP_
#define LONG_READ_MAPPER_HPP_

#include "pe_utils.hpp"
#include "graphio.hpp"
namespace path_extend {

class SimpleLongReadMapper {
public:
    SimpleLongReadMapper(conj_graph_pack& gp)
            : gp_(gp)/*,
              mapper_(gp_.g, gp_.index, gp_.kmer_mapper, gp_.k_value + 1),
              same_edge_corr_(gp_.g),
              ps_(gp_.g),
              gap_closer_(gp_.g, &ps_) */{
        //paths_searcher_config conf;
        //conf.depth_neigh_search = 5;  // max path len (in edges)
        //conf.max_len_path = 100000;  // max path len (in k-mers)
        //conf.max_num_vertices = 100;  // max number of visited vertices
        //ps_.Initialize(conf);
    }

    void ProcessSingleReadLibrary(
            const io::SequencingLibrary<debruijn_config::DataSetData>& lib,
            PathStorage<Graph>& storage) {
        if (cfg::get().use_multithreading) {
            auto single_streams = single_binary_readers(
                    lib, false, false);
            if (single_streams->size() == (size_t) 1) {
                ProcessReads(*single_streams, storage);
            } else {
                ProcessSingleReadsParallel(*single_streams, storage);
            }
        } else {
            auto single_streams = single_binary_readers(lib, false, false);
            single_streams->release();
            io::MultifileReader<io::SingleReadSeq> stream(single_streams->get(),
                                                          true);
            ProcessLib(stream, storage);
        }
    }

    template<class SingleRead>
    void ProcessLib(io::IReader<SingleRead>& stream,
                    PathStorage<Graph>& storage) {
        INFO("Processing single reads (takes a while)");
        while (!stream.eof()) {
            SingleRead r;
            stream >> r;
            vector<EdgeId> path = ProcessSingleRead(r);
            storage.AddPath(path, 1, true);
        }
    }
private:
    template<class SingleRead>
    void ProcessReads(
            io::ReadStreamVector<io::IReader<SingleRead> >& streams,
            PathStorage<Graph>& storage) {
        INFO("Mapping single reads (takes a while)");
        io::IReader<SingleRead>& stream = streams.back();
        stream.reset();
        while (!stream.eof()) {
            SingleRead r;
            stream >> r;
            vector<EdgeId> path = ProcessSingleRead(r);
            storage.AddPath(path, 1, true);
        }

    }

    template<class SingleRead>
    void ProcessSingleReadsParallel(
            io::ReadStreamVector<io::IReader<SingleRead> >& streams,
            PathStorage<Graph>& storage) {
        INFO("Mapping single reads (takes a while)");
        size_t nthreads = streams.size();
        vector<PathStorage<Graph>*> buffer_storages(nthreads);

        for (size_t i = 0; i < nthreads; ++i) {
            buffer_storages[i] = new PathStorage<Graph>(gp_.g);
        }

        size_t counter = 0;
        static const double coeff = 1.3;
        size_t count_mapped_reads = 0;
        size_t count_unmapped_reads = 0;
        size_t count_mapped_reads_size_one = 0;
        #pragma omp parallel num_threads(nthreads)
        {
            #pragma omp for reduction(+ : counter)
            for (size_t i = 0; i < nthreads; ++i) {
                size_t size = 0;
                size_t limit = 1000000;
                SingleRead r;
                io::IReader<SingleRead>& stream = streams[i];
                stream.reset();
                bool end_of_stream = false;
                DEBUG("Starting " << omp_get_thread_num());
                while (!end_of_stream) {
                    end_of_stream = stream.eof();

                    while (!end_of_stream && size < limit) {
                        stream >> r;
                        ++counter;
                        ++size;
                        vector<EdgeId> path = ProcessSingleRead(r);
                        buffer_storages[i]->AddPath(path, 1, true);
                        #pragma omp critical
                        {
                            if (path.size() == 0) {
                                count_unmapped_reads++;
                            } else {
                                count_mapped_reads++;
                                if (path.size() == 1) {
                                    count_mapped_reads_size_one++;
                                }
                            }
                        }
                        end_of_stream = stream.eof();
                    }

                    #pragma omp critical
                    {
                        DEBUG("Merging " << omp_get_thread_num() << " " << buffer_storages[i]->size());
                        storage.AddStorage(*(buffer_storages[i]));
                        DEBUG("New basket size " << storage.size());
                        DEBUG("Thread number " << omp_get_thread_num() << " is going to increase its limit by " << coeff << " times, current limit is " << limit);
                    }
                    buffer_storages[i]->Clear();
                    limit = coeff * limit;
                }
            }
            DEBUG("Thread number " << omp_get_thread_num() << " finished");
        }
        INFO("Count unmapped reads " << count_unmapped_reads
                          << " mapped reads " << count_mapped_reads
                          << " with size one " << count_mapped_reads_size_one);
        for (size_t i = 0; i < nthreads; ++i) {
            delete buffer_storages[i];
        }
        DEBUG("Done");
    }

    template<class SingleRead>
    vector<EdgeId> ProcessSingleRead(const SingleRead& r) {
        //TODO: if we can really use following code then we should delete everything about SameEdgeDeletionCorrector and CloseGapsCorrector!!!
        auto mapper = MapperInstance(gp_);
        return mapper->FindReadPath(r.sequence());
        /*MappingPath<EdgeId> path;
         path.join(mapper_.MapSequence(r.sequence()));
         SimpleMappingContig mc(r.sequence(), path);
         MappingContig * dc = same_edge_corr_.Correct(&mc);
         MappingContig * gc = gap_closer_.Correct(dc);
         return gc->PathSeq();*/
    }

    conj_graph_pack& gp_;
    //ExtendedSequenceMapper<Graph, conj_graph_pack::index_t> mapper_;
    //SameEdgeDeletionCorrector same_edge_corr_;
    //DijkstraSearcher ps_;
    //CloseGapsCorrector gap_closer_;
};

}/*path_extend*/


#endif /* LONG_READ_MAPPER_HPP_ */
