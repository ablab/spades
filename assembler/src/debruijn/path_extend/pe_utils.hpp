/*
 * pe_utils.hpp
 *
 *  Created on: Nov 27, 2012
 *      Author: andrey
 */

#ifndef PE_UTILS_HPP_
#define PE_UTILS_HPP_

#include "bidirectional_path.hpp"
#include "graph_pack.hpp"

using namespace debruijn_graph;

namespace path_extend {


class GraphCoverageMap: public PathListener {

public:
    typedef std::multiset <BidirectionalPath *> MapDataT;


protected:
    Graph& g_;

    std::map <EdgeId, MapDataT * > edgeCoverage_;

    MapDataT * empty_;

    virtual void EdgeAdded(EdgeId e, BidirectionalPath * path, int /*gap*/) {
        auto iter = edgeCoverage_.find(e);
        if (iter == edgeCoverage_.end()) {
            edgeCoverage_.insert(std::make_pair(e, new MapDataT()));
        }
        edgeCoverage_[e]->insert(path);
    }

    virtual void EdgeRemoved(EdgeId e, BidirectionalPath * path) {
        auto iter = edgeCoverage_.find(e);
        if (iter != edgeCoverage_.end()) {
            if (iter->second->count(path) == 0) {
                DEBUG("Error erasing path from coverage map");
            } else {
                auto entry = iter->second->find(path);
                iter->second->erase(entry);
            }
        }
    }

public:
    GraphCoverageMap(Graph& g) : g_(g), edgeCoverage_() {
        empty_ = new MapDataT();
    }

    GraphCoverageMap(Graph& g, PathContainer& paths) : g_(g), edgeCoverage_() {
        empty_ = new MapDataT();

        for (size_t i = 0; i < paths.size(); ++i) {
            for (size_t j = 0; j < paths.Get(i)->Size(); ++j) {
                EdgeAdded(paths.Get(i)->At(j), paths.Get(i), paths.Get(i)->GapAt(i));
            }
            for (size_t j = 0; j < paths.GetConjugate(i)->Size(); ++j) {
                EdgeAdded(paths.GetConjugate(i)->At(j), paths.GetConjugate(i), paths.GetConjugate(i)->GapAt(i));
            }
        }
    }

	void Subscribe(BidirectionalPath * path) {
		path->Subscribe(this);
		for (size_t i = 0; i < path->Size(); ++i) {
			BackEdgeAdded(path->At(i), path, path->GapAt(i));
		}
	}

    virtual void FrontEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        EdgeAdded(e, path, gap);
    }

    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        EdgeAdded(e, path, gap);
    }

    virtual void FrontEdgeRemoved(EdgeId e, BidirectionalPath * path) {
        EdgeRemoved(e, path);
    }

    virtual void BackEdgeRemoved(EdgeId e, BidirectionalPath * path) {
        EdgeRemoved(e, path);
    }

    MapDataT * GetEdgePaths(EdgeId e) const {
        auto iter = edgeCoverage_.find(e);
        if (iter != edgeCoverage_.end()) {
            return iter->second;
        }

        return empty_;
    }


    int GetCoverage(EdgeId e) const {
        return (int) GetEdgePaths(e)->size();
    }


    bool IsCovered(EdgeId e) const {
        return GetCoverage(e) > 0;
    }

    bool IsCovered(const BidirectionalPath& path) const {
        for (size_t i = 0; i < path.Size(); ++i) {
            if (!IsCovered(path[i])) {
                return false;
            }
        }
        return true;
    }

    int GetCoverage(const BidirectionalPath& path) const {
        if (path.Empty()) {
            return 0;
        }

        int cov = GetCoverage(path[0]);
        for (size_t i = 1; i < path.Size(); ++i) {
            int currentCov = GetCoverage(path[i]);
            if (cov > currentCov) {
                cov = currentCov;
            }
        }

        return cov;
    }

    std::set<BidirectionalPath*> GetCoveringPaths(EdgeId e) const {
        auto mapData = GetEdgePaths(e);
        return std::set<BidirectionalPath*>(mapData->begin(), mapData->end());

    }

    std::set<BidirectionalPath*> GetCoveringPaths(const BidirectionalPath& path) const {
        std::set<BidirectionalPath*> result;

        if (!path.Empty()) {
            MapDataT * data;
            data = GetEdgePaths(path.Front());

            result.insert(data->begin(), data->end());

            for (size_t i = 1; i < path.Size(); ++i) {
                data = GetEdgePaths(path[i]);

                std::set<BidirectionalPath*> dataSet;
                dataSet.insert(data->begin(), data->end());

                for (auto iter = result.begin(); iter != result.end(); ) {
                    auto next = iter;
                    ++next;
                    if (dataSet.count(*iter) == 0) {
                        result.erase(iter);
                    }
                    iter = next;
                }
            }
        }

        return result;
    }

    int GetUniqueCoverage(EdgeId e) const {
        return (int) GetCoveringPaths(e).size();
    }

    int GetUniqueCoverage(const BidirectionalPath& path) const {
        return (int) GetCoveringPaths(path).size();
    }

    std::map <EdgeId, MapDataT * >::const_iterator begin() const {
        return edgeCoverage_.begin();
    }

    std::map <EdgeId, MapDataT * >::const_iterator end() const {
        return edgeCoverage_.end();
    }

    // DEBUG

    void PrintUncovered() const {
        INFO("Uncovered edges");
        int s = 0;
        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (!IsCovered(*iter)) {
                INFO(g_.int_id(*iter) << " (" << g_.length(*iter) << ") ~ " << g_.int_id(g_.conjugate(*iter)) << " (" << g_.length(g_.conjugate(*iter)) << ")");
                s += 1;
            }
        }
        INFO("Uncovered edges " << s / 2);
    }

    void PrintMulticovered() const {
        INFO("Multicovered edges");
        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            auto paths = GetCoveringPaths(*iter);
            if (paths.size() > 1 && g_.length(*iter) > 1000) {
                INFO(g_.int_id(*iter) << " (" << g_.length(*iter) << "). " << " Covered: " << paths.size());
                for (auto path = paths.begin(); path != paths.end(); ++path) {
                    (*path)->Print();
                }
                INFO("=====");
            }
        }
    }
};



struct paths_searcher_config{
    size_t max_num_vertices;
    size_t depth_neigh_search;
    size_t max_len_path;
};

class PathsSearcher{

protected:
    Graph & g_;
    paths_searcher_config conf_;

public:
    PathsSearcher(Graph & g) : g_(g) {

    }

    virtual ~PathsSearcher() {

    }

    void Initialize(paths_searcher_config conf){
        conf_ = conf;
    }

    virtual map<VertexId, vector<EdgeId> > FindShortestPathsFrom(VertexId v) = 0;
};


class DijkstraSearcher : public PathsSearcher{

public:
    DijkstraSearcher(Graph & g) : PathsSearcher(g) {
    }

    map<VertexId, vector<EdgeId> > FindShortestPathsFrom(VertexId v){
        map<VertexId, vector<EdgeId> > short_paths;

        multimap<size_t, VertexId> dist_v;
        map<VertexId, size_t> v_dist;
        map<VertexId, size_t> v_depth;
        set<VertexId> visited;

        // insertion of the initial vertex
        vector<EdgeId> empty_path;
        dist_v.insert(pair<size_t, VertexId>(0, v));
        v_dist.insert(pair<VertexId, size_t>(v, 0));
        short_paths.insert(pair<VertexId, vector<EdgeId> >(v, empty_path));
        v_depth[v] = 0;

        size_t num_visited = 0;

        while((visited.size() < conf_.max_num_vertices) && (dist_v.size() != 0)) {

            VertexId cur_v = dist_v.begin()->second;
            size_t cur_dist = dist_v.begin()->first;

            size_t cur_depth;
            if(v_depth.find(cur_v) != v_depth.end()) {
                cur_depth = v_depth[cur_v];
            }
            else {
                size_t min_depth = 100000;
                bool is_defined = false;

                // defining of depth
                vector<EdgeId> in_edges = g_.IncomingEdges(cur_v);
                for(auto e = in_edges.begin(); e!= in_edges.end(); e++){
                    VertexId w = g_.EdgeStart(*e);
                    if(v_depth.find(w) != v_depth.end())
                        if(min_depth > v_depth[w]){
                            min_depth = v_depth[w];
                            is_defined = true;
                        }
                }

                if(is_defined){
                    cur_depth = min_depth + 1;
                }
                else{
                    cur_depth = 0;
                }
                v_depth[cur_v] = cur_depth;
            }

            if((cur_depth <= conf_.depth_neigh_search)){
                vector<EdgeId> out_edges = g_.OutgoingEdges(cur_v);

                for(auto e = out_edges.begin(); e != out_edges.end(); e++){
                    VertexId cur_neigh = g_.EdgeEnd(*e);

                    if(visited.find(cur_neigh) == visited.end()){
                        size_t new_neigh_dist = g_.length(*e) + cur_dist;
                        bool is_replaced = false;
                        if(v_dist.find(cur_neigh) != v_dist.end()){
                            size_t old_neigh_dist = v_dist[cur_neigh];

                            if(old_neigh_dist > new_neigh_dist){
                                is_replaced = true;

                                for(auto it = dist_v.find(old_neigh_dist); it != dist_v.end(); it++)
                                    if(it->second == cur_neigh){
                                        dist_v.erase(it);
                                        break;
                                    }
                            }
                        }
                        else {
                            is_replaced = true;
                        }

                        if(is_replaced && new_neigh_dist <= conf_.max_len_path){
                            dist_v.insert(pair<size_t, VertexId>(new_neigh_dist, cur_neigh));
                            v_dist[cur_neigh] = new_neigh_dist;

                            short_paths[cur_neigh] = short_paths[cur_v];
                            short_paths[cur_neigh].push_back(*e);
                        }
                    }
                }
            }
            else{
                break;
            }

            num_visited++;
            visited.insert(cur_v);

            // erasing of visited element;
            for(auto it = dist_v.find(cur_dist); it != dist_v.end(); it++){
                if(it->second == cur_v){
                    dist_v.erase(it);
                    v_dist.erase(it->second);
                    break;
                }
            }
        }

        return short_paths;
    }
};




class MappingContig{

public:
    virtual ~MappingContig() {

    }

    //virtual Sequence Seq() = 0;

    virtual vector<EdgeId> PathSeq() = 0;

    virtual size_t Size() = 0;

};


class SimpleMappingContig : public MappingContig{

    Sequence seq;

    MappingPath<EdgeId> map_path;

public:

    SimpleMappingContig(Sequence s, MappingPath<EdgeId> p) : seq(s), map_path(p) {}

    vector<EdgeId> PathSeq() {
        return map_path.simple_path().sequence();
    }

    size_t Size() {
        return map_path.size();
    }
};


class DeletedEdgeMappingContig : public MappingContig{

    MappingContig *c_;

    set<int> del_edges_;

    vector<EdgeId> path_seq;

    int new_size;

    void DefinePathSeq(){
        new_size = (int) c_->Size() - (int) del_edges_.size();

        VERIFY(new_size >= 0);

        vector<EdgeId> old_path = c_->PathSeq();
        for (size_t i = 0; i < old_path.size(); i++)
            if (del_edges_.find((int) i) == del_edges_.end())
                path_seq.push_back(old_path[i]);
    }

public:
    DeletedEdgeMappingContig(MappingContig * & c, set<int> del_edges) : c_(c), del_edges_(del_edges){
        new_size = -1;
    }

    DeletedEdgeMappingContig(MappingContig * & c, vector<int> del_edges) : c_(c){
        for (auto e = del_edges.begin(); e != del_edges.end(); e++)
            del_edges_.insert(*e);
        new_size = -1;
    }

    //Sequence Seq() { return c_->Seq(); }

    vector<EdgeId> PathSeq(){
        if (new_size == -1)
            DefinePathSeq();
        return path_seq;
    }

    size_t Size(){
        if (new_size == -1)
            DefinePathSeq();
        return new_size;
    }
};

class ReplacedPathMappingContig : public MappingContig{
    MappingContig *c_;

    vector<EdgeId> new_path_;

public:
    ReplacedPathMappingContig(MappingContig * c, vector<EdgeId> new_path) : c_(c), new_path_(new_path) { }
    //Sequence Seq() { return c_->Seq(); }
    vector<EdgeId> PathSeq() { return new_path_; }
    size_t Size() { return new_path_.size(); }
};


//----------------------------------------------------------------------------
class ContigStorage{
public:
    virtual void Add(MappingContig *new_contig) = 0;

    virtual size_t Size() = 0;

    virtual MappingContig* & operator[](int index) = 0;

    virtual void DeleteByIndex(int index) = 0;

    virtual void DeleteByIndexes(set<int> indexes) = 0;

    virtual void DeleteByIndexes(vector<int> indexes) = 0;

    virtual ~ContigStorage(){

    }
};

class SimpleContigStorage : public ContigStorage{
    vector<MappingContig *> stor;
    set<int> del_index_;
    vector<int> real_index;

    void RedefineRealIndexes(){
        for(size_t i = 0; i < stor.size(); i++)
            if(del_index_.find((int) i) == del_index_.end())
                real_index.push_back((int) i);
    }

public:
    void Add(MappingContig *new_contig) {
        stor.push_back(new_contig);
        real_index.push_back((int) real_index.size());
    }

    size_t Size() { return real_index.size(); }

    MappingContig* & operator[](int index){
        VERIFY(index >= 0 && index < (int) real_index.size());
        return stor[real_index[index]];
    }

    void DeleteByIndex(int index){
        real_index.clear();
        del_index_.insert(index);

        RedefineRealIndexes();
    }

    void DeleteByIndexes(set<int> indexes){
        real_index.clear();

        for(auto it = indexes.begin(); it != indexes.end(); it++)
            del_index_.insert(*it);
        RedefineRealIndexes();
    }

    void DeleteByIndexes(vector<int> indexes){
        real_index.clear();

        for(auto it = indexes.begin(); it != indexes.end(); it++)
            del_index_.insert(*it);
        RedefineRealIndexes();
    }

    virtual ~SimpleContigStorage(){
        for(auto it = stor.begin(); it != stor.end(); it++)
            delete *it;
    }
};



class ContigCorrector{
protected:
    Graph& g_;
public:
    ContigCorrector(Graph& g) : g_(g) {}

    virtual ~ContigCorrector() {

    }

    virtual ContigStorage * Correct(ContigStorage *storage) = 0;

    virtual MappingContig * Correct(MappingContig *contig) = 0;

};


class SameEdgeDeletionCorrector : public ContigCorrector{
public:
    SameEdgeDeletionCorrector(Graph &g) : ContigCorrector(g) {}

    ContigStorage * Correct(ContigStorage *contigs) {
        for (size_t i = 0; i < contigs->Size(); i++) {
            MappingContig *cur = (*contigs)[(int) i];
            (*contigs)[(int) i] = Correct(cur);
        }
        return contigs;
    }

    MappingContig * Correct(MappingContig *contig){
        vector<EdgeId> path = contig->PathSeq();
        set<int> red_edges;
        if (path.size() != 0) {
            EdgeId cur_edge = path[0];
            for (size_t i = 1; i < path.size(); i++) {
                EdgeId e = path[i];
                if (e != cur_edge){
                    cur_edge = e;
                }
                else {
                    red_edges.insert((int) i);
                }
            }
        }

        return new DeletedEdgeMappingContig(contig, red_edges);
    }
};


class CloseGapsCorrector : public ContigCorrector{

    PathsSearcher * ps_;
    vector<int> incorr_contigs;

    bool AreEdgesConnected(EdgeId e1, EdgeId e2){
        return g_.EdgeEnd(e1) == g_.EdgeStart(e2);
    }

    vector<EdgeId> ClosePathGap(vector<EdgeId> path, vector<size_t> gap_index){
        vector<EdgeId> new_path;
        size_t current_gap = 0;
        for(size_t i = 0; i < path.size() - 1; i++){

            EdgeId cur_edge = path[i];
            new_path.push_back(cur_edge);
            if(i == gap_index[current_gap]){
                VertexId start = g_.EdgeEnd(cur_edge);
                VertexId end = g_.EdgeStart(path[i + 1]);

                map<VertexId, vector<EdgeId> > path_map = ps_->FindShortestPathsFrom(start);

                if(path_map.find(end) != path_map.end()){
                    vector<EdgeId> add_path = path_map[end];

                    for(auto e = path_map[end].begin(); e != path_map[end].end(); e++)
                        if(g_.EdgeStart(*e) != g_.EdgeEnd(*e))
                            new_path.push_back(*e);
                }
                else
                    DEBUG("Gap from " + ToString(g_.int_id(start)) + " to " + ToString(g_.int_id(end))
                            + " won't closed");
                current_gap++;
            }
        }

        new_path.push_back(path[path.size() - 1]);
        return new_path;
    }

public:
    CloseGapsCorrector(Graph & g, PathsSearcher * ps) : ContigCorrector(g), ps_(ps) {
    }

    ContigStorage * Correct(ContigStorage *storage){
        for (size_t i = 0; i < storage->Size(); i++){
            auto cur_mc = (*storage)[(int) i];
            (*storage)[(int) i] = Correct(cur_mc);
        }
        return storage;
    }

    MappingContig * Correct(MappingContig *contig){
        vector<EdgeId> path = contig->PathSeq();

        if (path.size() <= 1)
            return contig;

        vector<size_t> gap_indexes;

        for (size_t i = 0; i < path.size() - 1; i++){
            EdgeId e1 = path[i];
            EdgeId e2 = path[i + 1];

            if (!AreEdgesConnected(e1, e2)){
                gap_indexes.push_back(i);
            }
        }

        if (gap_indexes.size() != 0){
            vector<EdgeId> new_path = ClosePathGap(path, gap_indexes);
            return new ReplacedPathMappingContig(contig, new_path);
        }
        else
            return contig;
    }
};


class DeleteUnconnectedContigCorrector : public ContigCorrector{

    bool IsPathConnected(vector<EdgeId> path){
        if(path.size() <= 1)
            return true;

        for(size_t i = 0; i < path.size() - 1; i++){
            EdgeId e1 = path[i], e2 = path[i + 1];
            if(g_.EdgeEnd(e1) != g_.EdgeStart(e2))
                return false;
        }

        return true;
    }

public:
    DeleteUnconnectedContigCorrector(Graph& g) : ContigCorrector(g) {}

    ContigStorage * Correct(ContigStorage *storage){

        vector<int> incorr_contigs;
        for(int i = 0; i < (int) storage->Size(); i++){
            auto path = (*storage)[i]->PathSeq();
            if(!IsPathConnected(path))
                incorr_contigs.push_back(i);
        }

        storage->DeleteByIndexes(incorr_contigs);
        return storage;
    }

    MappingContig * Correct(MappingContig *contig){
        return contig;
    }
};



class CompositeCorrector : public ContigCorrector{

    vector<ContigCorrector*> coll;
public:
    CompositeCorrector(Graph &g) : ContigCorrector(g) {}

    void AddCorrector(ContigCorrector* new_corr){
        coll.push_back(new_corr);
    }

    ContigStorage * Correct(ContigStorage *storage){
        ContigStorage *cur_stor = storage;
        for(auto it = coll.begin(); it != coll.end(); it++){
            cur_stor = (*it)->Correct(cur_stor);
        }
        return cur_stor;
    }

    MappingContig * Correct(MappingContig *contig){
        MappingContig *cur_contig = contig;
        for(auto it = coll.begin(); it != coll.end(); it++){
            cur_contig = (*it)->Correct(cur_contig);
        }
        return cur_contig;
    }

};


vector<Sequence> SplitContig(const io::SingleRead& read) {
    std::vector<Sequence> buffer;

    for(size_t i = 0; i < read.size(); i++) {
        size_t j = i;
        while(j < read.size() && is_nucl(read.GetSequenceString()[j])) {
            j++;
        }
        if(j > i) {
            buffer.push_back(Sequence(read.Substr(i, j)));
            i = j - 1;
        }
    }
    return buffer;
}

vector< vector<Sequence> > ReadContigsFromFile(const string& filename) {
    io::Reader reader(filename);
    io::SingleRead read;
    vector< vector<Sequence> > result;

    while (!reader.eof()) {
        reader >> read;
        result.push_back(SplitContig(read));
    }
    return result;
}


PathContainer CreatePathsFromContigs(conj_graph_pack & gp, const string& filename){
    DEBUG("Reading contigs");
    vector< vector<Sequence> > contigs = ReadContigsFromFile(filename);
    INFO(ToString(contigs.size()) + " contigs were read");

    vector< MappingPath<EdgeId> > paths;
    auto mapper = MapperInstance(gp);

    INFO("-- Paths before simplification");

    for(auto it = contigs.begin(); it != contigs.end(); it++) {
        MappingPath<EdgeId> path;
        for (auto seq = it->begin(); seq != it->end(); ++seq) {
            path.join(mapper->MapSequence(*seq));
        }
        if (path.size() == 0) {
            DEBUG("Empty path from");
            for (auto seq = it->begin(); seq != it->end(); ++seq) {
                DEBUG(seq->str());
            }
        }

        paths.push_back(path);
        for(size_t i = 0; i < path.size(); i++){
        	cout << gp.g.int_id(path[i].first) << " ";
        }
        cout << endl;
    }


    ContigStorage * storage = new SimpleContigStorage();
    for(size_t i = 0; i < contigs.size(); i++){
        storage->Add(new SimpleMappingContig(contigs[i][0], paths[i]));
    }

    DEBUG("Deleting repeat edges");
    SameEdgeDeletionCorrector * sedCorrector = new SameEdgeDeletionCorrector(gp.g);
    sedCorrector->Correct(storage);

    DEBUG("Closing gaps");
    PathsSearcher * ps = new DijkstraSearcher(gp.g);
    paths_searcher_config conf;
    conf.depth_neigh_search = 5; // max path len (in edges)
    conf.max_len_path = 100000;  // max path len (in k-mers)
    conf.max_num_vertices = 100; // max number of visited vertices
    ps->Initialize(conf);

    CloseGapsCorrector * cgCorrector = new CloseGapsCorrector(gp.g, ps);
    cgCorrector->Correct(storage);


    DEBUG("Making paths from " << storage->Size() << " contigs");
    PathContainer bidirectionalPaths;
    for (size_t i = 0; i < storage->Size(); ++i ) {
        if ((*storage)[(int) i]->PathSeq().size() == 0)
            continue;

        BidirectionalPath * p = new BidirectionalPath(gp.g, (*storage)[(int) i]->PathSeq());
        BidirectionalPath * cp = new BidirectionalPath(p->Conjugate());
        bidirectionalPaths.AddPair(p, cp);
    }

    INFO("Created " << bidirectionalPaths.size() << " paths");
    return bidirectionalPaths;
}


class ScaffoldBreaker {

private:

    int min_gap_;

    PathContainer container_;

    void SplitPath(const BidirectionalPath& path) {
        size_t i = 0;

        while (i < path.Size()) {
            BidirectionalPath * p = new BidirectionalPath(path.graph(), path[i]);
            size_t rc_id = path.graph().int_id(path.graph().conjugate(path[i]));
            ++i;

            while(i < path.Size() and path.GapAt(i) <= min_gap_) {
                p->PushBack(path[i], path.GapAt(i));
                ++i;
            }

            BidirectionalPath * cp = new BidirectionalPath(p->Conjugate());
            cp->SetId(rc_id);
            container_.AddPair(p, cp);
        }
    }

public:

    ScaffoldBreaker(int min_gap): min_gap_(min_gap) {

    }

    void Split(PathContainer& paths) {
        for (auto it = paths.begin(); it != paths.end(); ++it) {
            SplitPath(*it.get());
        }
    }


    void clear() {
        container_.clear();
    }

    PathContainer& container() {
        return container_;
    }

};

}

#endif /* PE_UTILS_HPP_ */
