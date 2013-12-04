//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * pe_io.hpp
 *
 *  Created on: Mar 15, 2012
 *      Author: andrey
 */

#ifndef PE_IO_HPP_
#define PE_IO_HPP_


#include "bidirectional_path.hpp"
#include "io/osequencestream.hpp"

namespace path_extend {

using namespace debruijn_graph;

class ContigWriter {
protected:
    DECL_LOGGER("PathExtendIO")

protected:
	const Graph& g_;

    size_t k_;

	string ToString(const BidirectionalPath& path) const{
		stringstream ss;
		if (!path.Empty()) {
			ss <<g_.EdgeNucls(path[0]).Subseq(0, k_).str();
		}

		for (size_t i = 0; i < path.Size(); ++i) {
			int gap = i == 0 ? 0 : path.GapAt(i);
			if (gap > (int) k_) {
				for (size_t j = 0; j < gap - k_; ++j) {
					ss << "N";
				}
				ss << g_.EdgeNucls(path[i]).str();
			} else {
				int overlapLen = (int) k_ - gap;
				if (overlapLen >= (int) g_.length(path[i]) + (int) k_) {
					continue;
				}

				ss << g_.EdgeNucls(path[i]).Subseq(overlapLen).str();
			}
		}
		return ss.str();
	}

    Sequence ToSequence(const BidirectionalPath& path) const {
        SequenceBuilder result;

        if (!path.Empty()) {
            result.append(g_.EdgeNucls(path[0]).Subseq(0, k_));
        }
        for (size_t i = 0; i < path.Size(); ++i) {
            result.append(g_.EdgeNucls(path[i]).Subseq(k_));
        }

        return result.BuildSequence();
    }

    void FindPathsOrder(PathContainer& paths, multimap<VertexId, BidirectionalPath*>& starting) const {
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            BidirectionalPath* path = iter.get();
            starting.insert(make_pair(g_.EdgeStart(path->Front()), path));

            path = iter.getConjugate();
            starting.insert(make_pair(g_.EdgeStart(path->Front()), path));
        }
    }

    void ConstructFASTG(PathContainer& paths,
            map<BidirectionalPath*, string >& ids,
            map<BidirectionalPath*, vector<string> >& next_ids) {

        int counter = 1;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            BidirectionalPath* p = iter.get();
            string name = io::MakeContigId(counter++, p->Length() + k_, p->Coverage(), p->GetId());
            ids.insert(make_pair(p, name));
            ids.insert(make_pair(iter.getConjugate(), name + "'"));
        }

        multimap<VertexId, BidirectionalPath*> starting;
        set<VertexId> visited;
        queue<BidirectionalPath*> path_queue;
        FindPathsOrder(paths, starting);

        auto it = starting.begin();
        for (; it != starting.end(); ++it) {
            VertexId v = it->first;
            if (visited.count(v) > 0)
                continue;

            visited.insert(v);
            visited.insert(g_.conjugate(v));
            while (it != starting.upper_bound(v)) {
                path_queue.push(it->second);
                ++it;
            }

            while (!path_queue.empty()) {
                BidirectionalPath* p = path_queue.front();
                path_queue.pop();
                next_ids.insert(make_pair(p, vector<string>()));

                v = g_.EdgeEnd(p->Back());

                bool add_new_paths = (visited.count(v) == 0);
                for (auto v_it = starting.lower_bound(v); v_it != starting.upper_bound(v); ++v_it) {
                    BidirectionalPath* next_path = v_it->second;
                    next_ids[p].push_back(ids[p]);
                    TRACE("Node " << counter);
                    if (add_new_paths)
                        path_queue.push(next_path);
                }
                if (add_new_paths) {
                    visited.insert(v);
                    visited.insert(g_.conjugate(v));
                }
            }
        }
    }


public:
    ContigWriter(const Graph& g): g_(g), k_(g.k()){

    }

    void writeEdges(const string& filename) const {
        INFO("Outputting edges to " << filename);
        io::osequencestream_with_data_for_scaffold oss(filename);

        set<EdgeId> included;
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (included.count(*iter) == 0) {
                oss.setCoverage(g_.coverage(*iter));
                oss.setID((int) g_.int_id(*iter));
                oss << g_.EdgeNucls(*iter);

                included.insert(*iter);
                included.insert(g_.conjugate(*iter));
            }
        }
        DEBUG("Contigs written");
    }


    void writePathEdges(PathContainer& paths, const string& filename) const {
		INFO("Outputting path data to " << filename);
		std::ofstream oss;
        oss.open(filename.c_str());
        int i = 1;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
			oss << i << endl;
			i++;
            BidirectionalPath* path = iter.get();
            if (path->GetId() % 2 != 0) {
                path = path->GetConjPath();
            }
            oss << "PATH " << path->GetId() << " " << path->Size() << " " << path->Length() + k_ << endl;
            for (size_t j = 0; j < path->Size(); ++j) {
			    oss << g_.int_id(path->At(j)) << " " << g_.length(path->At(j)) << endl;
            }
            oss << endl;
		}
		oss.close();
		DEBUG("Edges written");
	}

    void writePaths(PathContainer& paths, const string& filename) const {

        INFO("Writing contigs to " << filename);
        io::osequencestream_with_data_for_scaffold oss(filename);
        int i = 0;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        	if (iter.get()->Length() <= 0){
        		continue;
        	}
        	DEBUG("NODE " << ++i);
            BidirectionalPath* path = iter.get();
            if (path->GetId() % 2 != 0) {
                path = path->GetConjPath();
            }
            path->Print();
        	oss.setID((int) path->GetId());
            oss.setCoverage(path->Coverage());
            oss << ToString(*path);
        }
        DEBUG("Contigs written");
    }

    void WritePathsToFASTG(PathContainer& paths, const string& filename, const string& fastafilename) {
        map<BidirectionalPath*, string > ids;
        map<BidirectionalPath*, vector<string> > next_ids;

        INFO("Constructing FASTG file from paths " << filename);
        ConstructFASTG(paths, ids, next_ids);

        INFO("Writing contigs in FASTG to " << filename);
        INFO("Writing contigs in FASTQ to " << fastafilename);
        io::osequencestream_for_fastg fastg_oss(filename);
        io::osequencestream_with_id oss(fastafilename);
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            BidirectionalPath* path = iter.get();
            if (path->Length() <= 0){
                continue;
            }
            DEBUG(ids[path]);

            oss.setID((int) path->GetId());
            oss.setCoverage(path->Coverage());
            oss << ToString(*path);
            fastg_oss.set_header(ids[path]);
            fastg_oss << next_ids[path] << ToString(*path);
        }
    }
};


class PathInfoWriter {
protected:
    DECL_LOGGER("PathExtendIO")


public:

    void writePaths(PathContainer& paths, const string& filename){
        std::ofstream oss(filename.c_str());

        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            iter.get()->Print(oss);
        }

        oss.close();
    }
};

}

#endif /* PE_IO_HPP_ */
