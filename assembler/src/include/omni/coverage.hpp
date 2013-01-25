//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * coverage.hpp
 *
 *  Created on: Jun 21, 2011
 *      Author: sergey
 */

#ifndef COVERAGE_HPP_
#define COVERAGE_HPP_

#include <unordered_map>
#include "logger/logger.hpp"
#include "io/reader.hpp"
#include "io/read_stream_vector.hpp"
#include "perfcounter.hpp"

namespace omnigraph {

template<class Graph>
class CoverageIndex: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef unordered_map<EdgeId, int> map_type;

private:

//	map_type storage_;

	size_t KPlusOneMerCoverage(EdgeId edge) const {
		return (size_t) (coverage(edge) * this->g().length(edge));
	}


	template<class ReadThreader>
	Path<EdgeId> ProcessSequence(const ReadThreader& threader, const Sequence& sequence) const {
        return threader.MapSequence(sequence);
    }

    void AddPathsToGraph(const Path<EdgeId>& path) {

        if (path.sequence().size() == 0) 
            return;

        const vector<EdgeId>& edges_list = path.sequence();

        for (auto it = edges_list.cbegin(); it != edges_list.cend(); ++it) {
            IncCoverage(*it, this->g().length(*it));
        }
        IncCoverage(edges_list[0], -int(path.start_pos()));
        EdgeId last = edges_list[edges_list.size() - 1];
        IncCoverage(last, int(path.end_pos()) - int(this->g().length(last)));
	}

    void IncCoverageInMap(EdgeId edge, int toAdd, map_type& map) {
        //VERIFY(toAdd >= 0);
        map[edge] += toAdd;
        VERIFY(map[edge] >= 0);
    }

    void AddPathsToMap(const Path<EdgeId>& path, map_type& map) {

        if (path.sequence().size() == 0)
            return;

        const vector<EdgeId>& edges_list = path.sequence();

        for (auto it = edges_list.cbegin(); it != edges_list.cend(); ++it) {
            IncCoverageInMap(*it, this->g().length(*it), map);
        }
        IncCoverageInMap(edges_list[0], -int(path.start_pos()), map);
        EdgeId last = edges_list[edges_list.size() - 1];
        IncCoverageInMap(last, int(path.end_pos()) - int(this->g().length(last)), map);
    }

public:
	CoverageIndex(Graph &g) :
		GraphActionHandler<Graph> (g, "CoverageIndex") {
	}

	virtual ~CoverageIndex() {
	}

	void SetCoverage(EdgeId edge, int cov) {
		VERIFY(cov >= 0);
		edge->SetCoverage(cov);
	}

	/**
	 * Returns average coverage of the edge
	 */
	double coverage(EdgeId edge) const {
		return (double) edge->GetRawCoverage() / this->g().length(edge);
	}

	/**
	 * Method increases coverage value
	 */
	void IncCoverage(EdgeId edge, int toAdd) {
		edge->IncCoverage(toAdd);
		VERIFY(edge->GetRawCoverage() >= 0);
	}

	/**
	 * Method increases coverage value by 1
	 */
	void IncCoverage(EdgeId edge) {
		IncCoverage(edge, 1);
	}

    template<class ReadThreader, class Read>
    void FillIndex(io::IReader<Read>& stream, const ReadThreader& threader) {

        INFO("Processing reads (takes a while)");
        size_t counter = 0;
        stream.reset();

        while (!stream.eof()) {
            Read r;
            stream >> r;
            Path<EdgeId> path = ProcessSequence(threader, r.sequence());
            AddPathsToGraph(path);

            VERBOSE_POWER(++counter, " reads processed");
        }

        INFO("DeBruijn graph coverage counted, reads used: " << counter);
    }

	template<class ReadThreader, class Read>
	void FillParallelIndex(io::ReadStreamVector< io::IReader<Read> >& streams, const ReadThreader& threader, size_t buffer_size) {

        INFO("Processing reads (takes a while)");
        perf_counter pc;
        size_t counter = 0;

        size_t nthreads = streams.size();
        size_t buf_size = buffer_size / (nthreads * (sizeof(Path<EdgeId>) + 32) );

        #pragma omp parallel num_threads(nthreads)
        {
            #pragma omp for reduction(+ : counter)
            for (size_t i = 0; i < nthreads; ++i) {

                Read r;
                io::IReader<Read>& stream = streams[i];
                stream.reset();
                std::vector< Path<EdgeId> > buffer(buf_size);

                size_t j = 0;
                while (!stream.eof()) {
                    stream >> r;
                    ++counter;
                    buffer[j++] = ProcessSequence(threader, r.sequence());

                    if (j == buf_size) {
                        j = 0;

                        #pragma omp critical
                        {
                            for (size_t l = 0; l < buf_size; ++l) {
                                AddPathsToGraph(buffer[l]);
                            }
                        }
                    }
                }

                #pragma omp critical
                {
                    for (size_t l = 0; l < j; ++l) {
                        AddPathsToGraph(buffer[l]);
                    }
                }
            }

        }

        INFO("DeBruijn graph coverage counted, reads used: " << counter);

        INFO("Elapsed time: " << pc.time_ms());
	}

    template<class ReadThreader, class Read>
    void FillFastParallelIndex(io::ReadStreamVector< io::IReader<Read> >& streams, const ReadThreader& threader) {

        INFO("Processing reads (takes a while)");
        perf_counter pc;
        size_t counter = 0;

        size_t nthreads = streams.size();
//
        std::vector<map_type*> maps(nthreads);
//        maps[0] = &storage_;

        for (size_t i = 0; i < nthreads; ++i) {
            maps[i] = new map_type();
        }


        #pragma omp parallel num_threads(nthreads)
        {
            #pragma omp for reduction(+ : counter)
            for (size_t i = 0; i < nthreads; ++i) {

                Read r;
                io::IReader<Read>& stream = streams[i];
                stream.reset();
                Path<EdgeId> path;

                while (!stream.eof()) {
                    stream >> r;
                    ++counter;
                    path = ProcessSequence(threader, r.sequence());

                    AddPathsToMap(path, *maps[i]);
                }
            }
        }

        INFO("Merging maps");
        for (size_t i = 0; i < nthreads; ++i) {
            for (auto it = maps[i]->begin(); it != maps[i]->end(); ++it) {
            	it->first->IncCoverage(it->second);
            }
            delete maps[i];
        }

        INFO("DeBruijn graph coverage counted, reads used: " << counter);

        INFO("Elapsed time: " << pc.time_ms());
    }


	virtual void HandleDelete(EdgeId edge) {
		SetCoverage(edge, 0);
	}

	virtual void HandleMerge(const vector<EdgeId>& oldEdges, EdgeId newEdge) {
		size_t coverage = 0;
		for (auto it = oldEdges.begin(); it
				!= oldEdges.end(); ++it) {
			coverage += KPlusOneMerCoverage(*it);
		}
		SetCoverage(newEdge, coverage);
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		IncCoverage(new_edge, KPlusOneMerCoverage(edge2));
		IncCoverage(new_edge, KPlusOneMerCoverage(edge1));
	}

	virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
//		size_t length1 = this->g().length(newEdge1);
//		size_t length = this->g().length(oldEdge);
//		size_t coverage = KPlusOneMerCoverage(oldEdge);
//		size_t coverage1 = coverage * length1 / length;
//		if (coverage1 == 0)
//			coverage1 = 1;
//		size_t coverage2 = coverage - coverage1;
//		if (coverage2 == 0)
//			coverage2 = 1;
//		SetCoverage(newEdge1, coverage1);
//		SetCoverage(newEdge2, coverage2);
		double avg_cov = coverage(oldEdge);
		SetCoverage(newEdge1, size_t(max(1., math::round(avg_cov * this->g().length(newEdge1)))));
		SetCoverage(newEdge2, size_t(max(1., math::round(avg_cov * this->g().length(newEdge2)))));
	}

 	void HandleVertexSplit(VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges, vector<double> &split_coefficients, VertexId oldVertex) {
		 DEBUG("HandleMerge by coverage handler");
 		 size_t n = newEdges.size();
		 for(size_t j = 0; j < n; j++) {
			 EdgeId old_ID = newEdges[j].first;
			 EdgeId new_ID = newEdges[j].second;
			 IncCoverage(new_ID, floor(KPlusOneMerCoverage(old_ID)*split_coefficients[j]));
		 }
 	 }

 	/*
 	 * Is thread safe if edges different threads process different edges.
 	 */
	bool IsThreadSafe() const {
		return true;
	}
};

}

#endif /* COVERAGE_HPP_ */
