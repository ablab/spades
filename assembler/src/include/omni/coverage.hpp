//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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

#include <tr1/unordered_map>
#include "logger/logger.hpp"
#include "io/reader.hpp"
#include "perfcounter.hpp"


namespace omnigraph {

template<class Graph>
class CoverageIndex: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef map<EdgeId, int> map_type;

private:

	map_type storage_;

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

public:
	CoverageIndex(Graph &g) :
		GraphActionHandler<Graph> (g, "CoverageIndex") {
	}

	virtual ~CoverageIndex() {
	}

	void SetCoverage(EdgeId edge, int cov) {

		VERIFY(cov >= 0);

		storage_[edge] = cov;

		VERIFY(storage_[edge] >= 0);
	}

	/**
	 * Returns average coverage of the edge
	 */
	double coverage(EdgeId edge) const {
		auto it = storage_.find(edge);
		if (it == storage_.end()) {
			return 0;
		}
		return (double) it->second / this->g().length(edge);
	}

	/**
	 * Method increases coverage value
	 */
	void IncCoverage(EdgeId edge, int toAdd) {
		//VERIFY(toAdd >= 0);
		storage_[edge] += toAdd;
		VERIFY(storage_[edge] >= 0);
	}

	/**
	 * Method increases coverage value by 1
	 */
	void IncCoverage(EdgeId edge) {
		IncCoverage(edge, 1);
	}

	template<class ReadThreader>
	void FillIndex(io::IReader<io::SingleRead>& reads_stream, const ReadThreader& threader) {
        
        typedef Path<EdgeId> DataType;
		
        INFO("Processing reads (takes a while)");
        
        perf_counter pc;
		
        size_t counter = 0;

        size_t buf_size = 4000000;
        size_t nthreads = 4;
        cout << "# threads" << endl;
        cin >>  nthreads;

        cout << "buf size" << endl;
        cin >> buf_size;

        //initialization
        vector<vector<Sequence> > vector_seq(nthreads, vector<Sequence>());

        for (size_t i = 0; i < nthreads; ++i) 
            vector_seq[i].reserve(buf_size);

        vector<vector<DataType>> buffer(nthreads, vector<DataType>());

        for (size_t i = 0; i < nthreads; ++i) 
            buffer[i].reserve(buf_size);

        std::vector<size_t> sizes(nthreads, 0);

		io::SingleRead read;
     

        //first iteration fill vector with push_back
        if (!reads_stream.eof()) {

            INFO("Filling Sequence Buffer");

            for (size_t j = 0; j < buf_size; ++j) {
                for(size_t i = 0; i < nthreads && !reads_stream.eof(); ++i) {
                    reads_stream >> read;
                    vector_seq[i].push_back(read.sequence());
                    sizes[i]++;
                    VERBOSE_POWER(++counter, " reads processed");
                }
            }

            INFO("Finished Filling Sequence Buffer");

            INFO("Filling Path Buffer");

            #pragma omp parallel num_threads(nthreads)
            {

            #pragma omp for 
            for (size_t i = 0; i < nthreads; ++i) {
                for (size_t j = 0; j < sizes[i]; ++j) {
                    Path<EdgeId> path = ProcessSequence(threader, vector_seq[i][j]);
                    buffer[i].push_back(path); 
                }
            }
            
            }
            INFO("Finished Filling Path Buffer");

            INFO("Threading Paths");

            for (size_t i = 0; i < nthreads; ++i) {
                for (size_t j = 0; j < sizes[i]; ++j) {
                    AddPathsToGraph(buffer[i][j]); 
                }
            }

            INFO("Finished Threading Paths");
            std::fill(sizes.begin(), sizes.end(), 0);
        }

        // now we are filling with operator []
        while (!reads_stream.eof()) {
            INFO("Filling Sequence Buffer");

            for (size_t j = 0; j < buf_size; ++j) {
                for(size_t i = 0; i < nthreads && !reads_stream.eof(); ++i) {
                    reads_stream >> read;
                    vector_seq[i][j] = read.sequence();
                    sizes[i]++;
                    VERBOSE_POWER(++counter, " reads processed");
                }
            }

            INFO("Finished Filling Sequence Buffer");

            INFO("Filling Path Buffer");

            #pragma omp parallel num_threads(nthreads)
            {

            #pragma omp for 
            for (size_t i = 0; i < nthreads; ++i) {
                for (size_t j = 0; j < sizes[i]; ++j) {
                    Path<EdgeId> path = ProcessSequence(threader, vector_seq[i][j]);
                    buffer[i][j] = path; 
                }
            }
            
            }
            INFO("Finished Filling Path Buffer");

            INFO("Threading Paths");

            for (size_t i = 0; i < nthreads; ++i) {
                for (size_t j = 0; j < buffer[i].size(); ++j) {
                    AddPathsToGraph(buffer[i][j]); 
                }
            }

            INFO("Finished Threading Paths");
            std::fill(sizes.begin(), sizes.end(), 0);
		}

        INFO("Elapsed time: " << pc.time_ms());
	}

	virtual void HandleDelete(EdgeId edge) {
		storage_.erase(edge);
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
};

}

#endif /* COVERAGE_HPP_ */
