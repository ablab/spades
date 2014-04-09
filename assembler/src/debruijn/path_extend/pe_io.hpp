//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
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

	string ToString(const BidirectionalPath& path) const {
		stringstream ss;
		if (path.IsInterstrandBulge() && path.Size() == 1) {
		    ss << g_.EdgeNucls(path.Back()).Subseq(k_, g_.length(path.Back())).str();
		    return ss.str();
		}

		if (!path.Empty()) {
			ss << g_.EdgeNucls(path[0]).Subseq(0, k_).str();
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

    void MakeIDS(PathContainer& paths,
                 map<BidirectionalPath*, string >& ids,
                 map<BidirectionalPath*, vector<string> >& next_ids) const {
        int counter = 1;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (iter.get()->Size() == 0)
                continue;

            BidirectionalPath* p = iter.get();
            BidirectionalPath* cp = iter.getConjugate();
            string name = io::MakeContigId(counter++, p->Length() + k_, p->Coverage(), p->GetId());
            ids.insert(make_pair(p, name));
            ids.insert(make_pair(cp, name + "'"));
            next_ids.insert(make_pair(p, vector<string>()));
            next_ids.insert(make_pair(cp, vector<string>()));
        }
    }

    void FindPathsOrder(PathContainer& paths,
                        multimap<VertexId, BidirectionalPath*>& v_starting,
                        multimap<EdgeId, BidirectionalPath*>& e_starting) const {
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (iter.get()->Size() == 0)
                            continue;

            BidirectionalPath* path = iter.get();
            DEBUG(g_.int_id(g_.EdgeStart(path->Front())) << " -> " << path->Size() << ", " << path->Length());
            v_starting.insert(make_pair(g_.EdgeStart(path->Front()), path));
            e_starting.insert(make_pair(path->Front(), path));


            path = iter.getConjugate();
            DEBUG(g_.int_id(g_.EdgeStart(path->Front())) << " -> " << path->Size() << ", " << path->Length());
            v_starting.insert(make_pair(g_.EdgeStart(path->Front()), path));
            e_starting.insert(make_pair(path->Front(), path));
        }
    }

    void VerifyIDS(PathContainer& paths,
                 map<BidirectionalPath*, string >& ids,
                 map<BidirectionalPath*, vector<string> >& next_ids,
                 multimap<VertexId, BidirectionalPath*>& v_starting,
                 multimap<EdgeId, BidirectionalPath*>& e_starting) const {

        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (iter.get()->Size() == 0)
                continue;

            BidirectionalPath* p = iter.get();
            VertexId v = g_.EdgeEnd(p->Back());
            TRACE("Node " << ids[p] << " is followed by: ");
            for (auto v_it = v_starting.lower_bound(v); v_it != v_starting.upper_bound(v); ++v_it) {
                TRACE("Vertex: " << ids[v_it->second]);
                auto it = find(next_ids[p].begin(), next_ids[p].end(), ids[v_it->second]);
                VERIFY(it != next_ids[p].end());
            }

            EdgeId e = p->Back();
            TRACE("Node " << ids[p] << " is followed by: ");
            for (auto e_it = e_starting.lower_bound(e); e_it != e_starting.upper_bound(e); ++e_it) {
               TRACE("Edge: " << ids[e_it->second]);
               auto it = find(next_ids[p].begin(), next_ids[p].end(), ids[e_it->second]);
               VERIFY(it != next_ids[p].end());
            }

            p = iter.getConjugate();
            v = g_.EdgeEnd(p->Back());
            TRACE("Node " << ids[p] << " is followed by: ");
            for (auto v_it = v_starting.lower_bound(v); v_it != v_starting.upper_bound(v); ++v_it) {
                TRACE("Vertex: " << ids[v_it->second]);
                auto it = find(next_ids[p].begin(), next_ids[p].end(), ids[v_it->second]);
                VERIFY(it != next_ids[p].end());
            }

            e = p->Back();
            TRACE("Node " << ids[p] << " is followed by: ");
            for (auto e_it = e_starting.lower_bound(e); e_it != e_starting.upper_bound(e); ++e_it) {
               TRACE("Edge: " << ids[e_it->second]);
               auto it = find(next_ids[p].begin(), next_ids[p].end(), ids[e_it->second]);
               VERIFY(it != next_ids[p].end());
            }
        }
    }

    void ConstructFASTG(PathContainer& paths,
            map<BidirectionalPath*, string >& ids,
            map<BidirectionalPath*, vector<string> >& next_ids) const {

        MakeIDS(paths, ids, next_ids);

        multimap<VertexId, BidirectionalPath*> v_starting;
        multimap<EdgeId, BidirectionalPath*> e_starting;
        //set<VertexId> visited;
        //queue<BidirectionalPath*> path_queue;
        FindPathsOrder(paths, v_starting, e_starting);

        DEBUG("RESULT");
        for (auto it = v_starting.begin(); it != v_starting.end(); ++it){
            BidirectionalPath* path = it->second;
            DEBUG(g_.int_id(it->first) << " -> " << path->Size() << ", " << path->Length());
        }

        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (iter.get()->Size() == 0)
                            continue;

            BidirectionalPath* path = iter.get();
            VertexId v = g_.EdgeEnd(path->Back());
            TRACE("Node " << ids[path] << " is followed by: ");
            for (auto v_it = v_starting.lower_bound(v); v_it != v_starting.upper_bound(v); ++v_it) {
                BidirectionalPath* next_path = v_it->second;
                TRACE(ids[next_path]);
                next_ids[path].push_back(ids[next_path]);
            }

            EdgeId e = path->Back();
            TRACE("Node " << ids[path] << " is followed by: ");
            for (auto e_it = e_starting.lower_bound(e); e_it != e_starting.upper_bound(e); ++e_it) {
                BidirectionalPath* next_path = e_it->second;
                TRACE(ids[next_path]);
                next_ids[path].push_back(ids[next_path]);
            }

            path = iter.getConjugate();
            v = g_.EdgeEnd(path->Back());
            TRACE("Node " << ids[path] << " is followed by: ");
            for (auto v_it = v_starting.lower_bound(v); v_it != v_starting.upper_bound(v); ++v_it) {
                BidirectionalPath* next_path = v_it->second;
                TRACE(ids[next_path]);
                next_ids[path].push_back(ids[next_path]);
            }

            e = path->Back();
            TRACE("Node " << ids[path] << " is followed by: ");
            for (auto e_it = e_starting.lower_bound(e); e_it != e_starting.upper_bound(e); ++e_it) {
                BidirectionalPath* next_path = e_it->second;
                TRACE(ids[next_path]);
                next_ids[path].push_back(ids[next_path]);
            }
        }

        VerifyIDS(paths, ids, next_ids, v_starting, e_starting);
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
        oss << paths.size() << endl;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
			//oss << i << endl;
			i++;
            BidirectionalPath* path = iter.get();
            if (path->GetId() % 2 != 0) {
                path = path->GetConjPath();
            }
            oss << "PATH " << path->GetId() << " " << path->Size() << " " << path->Length() + k_ << endl;
            for (size_t j = 0; j < path->Size(); ++j) {
			    oss << g_.int_id(path->At(j)) << " " << g_.length(path->At(j)) <<  " " << path->GapAt(j) << endl;
            }
            //oss << endl;
		}
		oss.close();
		DEBUG("Edges written");
	}

    void loadPaths(PathContainer& paths,  GraphCoverageMap& cover_map, const string& filename) const {
        paths.clear();
        map<size_t, EdgeId> int_ids;
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            int_ids.insert(make_pair(g_.int_id(*iter), *iter));
        }

        std::ifstream iss;
        iss.open(filename);
        size_t psize;
        iss >> psize;
        for(size_t i = 0; i < psize && !iss.eof(); ++i) {
            string s;
            size_t id;
            size_t size;
            size_t len;
            iss >> s >> id >> size >> len;
            VERIFY(s == "PATH");

            BidirectionalPath * path = new BidirectionalPath(g_);
            BidirectionalPath * conjugatePath = new BidirectionalPath(g_);
            paths.AddPair(path, conjugatePath);
            path->Subscribe(&cover_map);
            conjugatePath->Subscribe(&cover_map);
            for (size_t j = 0; !iss.eof() && j < size; ++j) {
                size_t eid;
                size_t elen;
                int gap;
                iss >> eid >> elen >> gap;
                EdgeId edge = int_ids[eid];
                conjugatePath->PushBack(edge, gap);
                VERIFY(g_.length(edge) == elen);
            }
            VERIFY(path->Length() + k_ == len);
        }
        VERIFY(psize == paths.size());
        iss.close();
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

    void WritePathsToFASTG(PathContainer& paths, const string& filename, const string& fastafilename) const {
        map<BidirectionalPath*, string > ids;
        map<BidirectionalPath*, vector<string> > next_ids;

        INFO("Constructing FASTG file from paths ");
        ConstructFASTG(paths, ids, next_ids);

        INFO("Writing contigs in FASTG to " << filename);
        INFO("Writing contigs in FASTA to " << fastafilename);
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
            DEBUG("NODE " << ids[path]);
            path->Print();
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
