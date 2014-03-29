#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <set>

using namespace std;

using namespace debruijn_graph;

namespace dipspades {

class OverlapGraph{
	map<size_t, set<size_t> > in_edges_;
	map<size_t, set<size_t> > out_edges_;
	map<size_t, pair<bool, int> > label;

	set<pair<size_t,size_t> > edges_;

	map<pair<size_t,size_t>, size_t> weights_;

	set<size_t> vertices_;
	map<size_t, pair<size_t, size_t> > map_id_rc_id;
	map<size_t ,size_t> id_ind;

	void CheckAndRemoveIsolatedVertex(size_t v){
		if(IncomingVerticesCount(v) == 0 && OutgoingVerticesCount(v) == 0){
			vertices_.erase(v);
		}
	}

public:
	OverlapGraph(){}
	OverlapGraph(vector<size_t> vertices, vector<size_t> id, vector<size_t> rc_id,
			map<size_t, vector<size_t> > in_edges, map<size_t, vector<size_t> > in_weight,
			map<size_t, vector<size_t> > out_edges, map<size_t, vector<size_t> > out_weight) {

		InitializeVertexSet(vertices, id, rc_id);
		InitializeIncomingVertices(in_edges, in_weight);
		InitializeOutgoingVertices(out_edges, out_weight);
	}

	void Clear(){
		vertices_.clear();
		map_id_rc_id.clear();
		id_ind.clear();
		label.clear();

		in_edges_.clear();
		out_edges_.clear();
		weights_.clear();
		edges_.clear();
	}

	void InitializeVertexSet(vector<size_t> vertices, vector<size_t> id, vector<size_t> rc_id){

		VERIFY(vertices.size() == id.size());
		VERIFY(vertices.size() == rc_id.size());

		size_t size = vertices.size();
		for(size_t i = 0; i < size; i++){

			auto v = vertices[i];
			vertices_.insert(v);

			map_id_rc_id[v] = pair<size_t, size_t>(id[i], rc_id[i]);
			id_ind[id[i]] = v;

			label[v] = pair<bool, int>(false, -1);
		}
	}

	void InitializeIncomingVertices(map<size_t, vector<size_t> > in_edges,
			map<size_t, vector<size_t> > in_weight){

		VERIFY(in_edges.size() == in_weight.size());

		auto it_v = in_edges.begin();
		auto it_w = in_weight.begin();

		for(size_t i = 0; i < in_edges.size(); i++){

			auto v = it_v->first;
			auto w = it_w->first;

			VERIFY(v == w);

			auto v_vect = it_v->second;
			auto w_vect = it_w->second;

			VERIFY(v_vect.size() == w_vect.size());

			for(size_t j = 0; j < v_vect.size(); j++){
				AddNeighVertices(v_vect[j], v, w_vect[j]);
			}

			it_v++; it_w++;
		}
	}

	void InitializeOutgoingVertices(map<size_t, vector<size_t> > out_edges,
			map<size_t, vector<size_t> > out_weight){

		VERIFY(out_edges.size() == out_weight.size());

		auto it_v = out_edges.begin();
		auto it_w = out_weight.begin();

		for(size_t i = 0; i < out_edges.size(); i++){

			auto v = it_v->first;
			auto w = it_w->first;

			VERIFY(v == w);

			auto v_vect = it_v->second;
			auto w_vect = it_w->second;

			VERIFY(v_vect.size() == w_vect.size());

			for(size_t j = 0; j < v_vect.size(); j++){
				AddNeighVertices(v, v_vect[j], w_vect[j]);
			}
			it_v++; it_w++;
		}
	}

	void AddNeighVertices(size_t start, size_t end, size_t weight){

		if(vertices_.find(start) == vertices_.end())
			vertices_.insert(start);
		if(vertices_.find(end) == vertices_.end())
		vertices_.insert(end);

		if(edges_.find(make_pair(start, end)) == edges_.end())
			edges_.insert(make_pair(start, end));
		weights_[make_pair(start, end)] = weight;

		in_edges_[end].insert(start);
		out_edges_[start].insert(end);
	}

	vector<size_t> GetVerticesWithoutInEdges(){
		vector<size_t> res;
		for(auto v = vertices_.begin(); v != vertices_.end(); v++)
			if(in_edges_.find(*v) == in_edges_.end() || in_edges_[*v].size() == 0){
				res.push_back(*v);
				break;
			}
		return res;
	}

	set<size_t> IncomingVertices(size_t v){
		set<size_t> res;
		if(in_edges_.find(v) != in_edges_.end())
			res = in_edges_[v];
		return res;
	}

	size_t IncomingVerticesCount(size_t v){
		if(in_edges_.find(v) != in_edges_.end())
			return in_edges_[v].size();
		return 0;
	}

	set<size_t> OutgoingVertices(size_t v){
		set<size_t> res;
		if(out_edges_.find(v) != out_edges_.end())
			res = out_edges_[v];
		return res;
	}

	size_t OutgoingVerticesCount(size_t v){
		if(out_edges_.find(v) != out_edges_.end())
			return out_edges_[v].size();
		return 0;
	}

	pair<bool, int> GetLabel(size_t v){
		if(label.find(v) != label.end())
			return label[v];
		return make_pair(false, -1);
	}

	void SetLabel(size_t v, bool bool_label, int value){
		if(label.find(v) != label.end())
			label[v] = make_pair(bool_label, value);
	}

	size_t GetWeightOf(pair<size_t, size_t> edge){
		if(weights_.find(edge) != weights_.end())
			return weights_[edge];
		return 0;
	}

	size_t GetWeightOf(size_t start, size_t end){
		return GetWeightOf(pair<size_t, size_t>(start, end));
	}

	void DeleteVertex(size_t v){
		vertices_.erase(v);
		auto in_set = in_edges_[v];
		auto out_set = out_edges_[v];

		for(auto w = in_set.begin(); w != in_set.end(); w++){
			DeleteEdge(*w, v);
		}

		in_edges_.erase(v);

		for(auto w = out_set.begin(); w != out_set.end(); w++){
			DeleteEdge(v, *w);
		}

		in_edges_.erase(v);

		size_t id = IdByInd(v);
		map_id_rc_id.erase(v);
		id_ind.erase(id);
	}

	void DeleteEdge(size_t start, size_t end){
		if(out_edges_.find(start) != out_edges_.end())
			out_edges_[start].erase(end);
		if(in_edges_.find(end) != in_edges_.end())
			in_edges_[end].erase(start);
		weights_.erase(pair<size_t, size_t>(start, end));
		edges_.erase(pair<size_t, size_t>(start, end));

//		CheckAndRemoveIsolatedVertex(start);
//		CheckAndRemoveIsolatedVertex(end);
	}

	set<pair<size_t,size_t> > Edges(){
		auto res = edges_;
		return res;
	}

	set<size_t> Vertices(){
		auto res = vertices_;
		return res;
	}

	bool IsEdgeExist(size_t start, size_t end){
		if(edges_.find(pair<size_t, size_t>(start, end)) != edges_.end())
			return true;
		return false;
	}

	bool IsVertexExist(size_t v){
		return vertices_.find(v) != vertices_.end();
	}

	size_t VerticesCount(){
		return vertices_.size();
	}

	size_t EdgesCount(){
		return edges_.size();
	}

	size_t IdByInd(size_t v){
		if(map_id_rc_id.find(v) != map_id_rc_id.end())
			return map_id_rc_id[v].first;
		return size_t(-1);
	}

	size_t RCIdByInd(size_t v){
		if(map_id_rc_id.find(v) != map_id_rc_id.end())
			return map_id_rc_id[v].second;
		return size_t(-1);
	}

	size_t IndById(size_t id){
		if(id_ind.find(id) != id_ind.end())
			return id_ind[id];
		return size_t(-1);
	}

	size_t IndOfRC(size_t v){
		return IndById(RCIdByInd(v));
	}

	bool IsVertexValid(size_t v){
		return (v != size_t(-1));
	}

	bool IsVertexIsolated(size_t v){
		return OutgoingVerticesCount(v) == 0 && IncomingVerticesCount(v) == 0;
	}

	vector<size_t> GetIsolatedEdges(){
		vector<size_t> res;
		for(auto v = vertices_.begin(); v != vertices_.end(); v++)
			if(IsVertexIsolated(*v))
				res.push_back(*v);
		return res;
	}
};

//----------------------------------------------------------------------------------------------

class OGD_StartVerticesDefiner{
protected:
	OverlapGraph &g_;
public:
	OGD_StartVerticesDefiner(OverlapGraph &g) : g_(g) {
	}
	virtual vector<size_t> GetStartVertices() = 0;
	virtual ~OGD_StartVerticesDefiner() {
	}
};

class OGD_GetParametrizedStartvertex : public OGD_StartVerticesDefiner{
	size_t start_vertex_;
public:
	OGD_GetParametrizedStartvertex(OverlapGraph &g, size_t start_vertex) :
		OGD_StartVerticesDefiner(g), start_vertex_(start_vertex){

	}

	vector<size_t> GetStartVertices(){
		vector<size_t> res;
		res.push_back(start_vertex_);
		return res;
	}
};

class OGD_GetIsolatedAndVerticesWithoutIncoming : public OGD_StartVerticesDefiner{
public:
	OGD_GetIsolatedAndVerticesWithoutIncoming(OverlapGraph &g) : OGD_StartVerticesDefiner(g){

	}

	vector<size_t> GetStartVertices(){

//		cout << "OGD_GetIsolatedAndVerticesWithoutIncoming starts" << endl;
//		cout << g_.VerticesCount() << " vertices in OG" << endl;

		vector<size_t> res;

		if(g_.VerticesCount() == 0)
			return res;

		vector<size_t> isolated = g_.GetIsolatedEdges();
		vector<size_t> noincoming = g_.GetVerticesWithoutInEdges();

		if(isolated.size() != 0)
			for(auto v = isolated.begin(); v != isolated.end(); v++)
				res.push_back(*v);

		if(noincoming.size() != 0)
				res.push_back(*noincoming.begin());

		if(res.size() == 0){
			size_t any_vertex = *g_.Vertices().begin();
			res.push_back(any_vertex);
		}

		return res;
	}
};

//----------------------------------------------------------------------------------------------
class OGD_DirectionDefiner{
protected:
	OverlapGraph &g_;

public:
	OGD_DirectionDefiner(OverlapGraph &g) : g_(g){

	}

	virtual set<size_t> GetDirectedVertices(size_t vertex) = 0;
	virtual set<size_t> GetAntidirectedVertices(size_t vertex) = 0;
	virtual ~OGD_DirectionDefiner(){

	}
};

class OGD_OutgoingDirection : public OGD_DirectionDefiner{
public:
	OGD_OutgoingDirection(OverlapGraph &g) : OGD_DirectionDefiner(g){

	}

	set<size_t> GetDirectedVertices(size_t vertex){
		return g_.OutgoingVertices(vertex);
	}

	set<size_t> GetAntidirectedVertices(size_t vertex){
		return g_.IncomingVertices(vertex);
	}
};
//----------------------------------------------------------------------------------------------

class OGD_OneVertexProcesser{
protected:
	OverlapGraph &g_;
	OGD_DirectionDefiner *direction_definer_;
public:
	OGD_OneVertexProcesser(OverlapGraph &g, OGD_DirectionDefiner *direction_definer) : g_(g),
		direction_definer_(direction_definer){
	}
	virtual void ProcessVertex(size_t vertex, set<size_t> &visited, set<size_t> &queue,
			map<size_t, vector<size_t> > &paths) = 0;
	virtual ~OGD_OneVertexProcesser(){

	}
};

class OGD_SimpleProcessing : public OGD_OneVertexProcesser{
public:
	OGD_SimpleProcessing(OverlapGraph &g, OGD_DirectionDefiner *direction_definer) :
		OGD_OneVertexProcesser(g, direction_definer){

	}

	void ProcessVertex(size_t vertex, set<size_t> &visited, set<size_t> &queue,
				map<size_t, vector<size_t> > &paths){
		if(visited.find(vertex) != visited.end())
			return;

		visited.insert(vertex);
		queue.erase(vertex);

		auto vert_for_visit = direction_definer_->GetDirectedVertices(vertex);
		for(auto neigh = vert_for_visit.begin(); neigh != vert_for_visit.end(); neigh++)
			if(visited.find(*neigh) == visited.end()){
				paths[*neigh] = paths[vertex];
				paths[*neigh].push_back(*neigh);

				queue.insert(*neigh);
			}
	}
};

class OGD_UniquePathProcessing : public OGD_OneVertexProcesser{
public:
	OGD_UniquePathProcessing(OverlapGraph &g, OGD_DirectionDefiner *direction_definer) :
		OGD_OneVertexProcesser(g, direction_definer){

	}

	void ProcessVertex(size_t vertex, set<size_t> &visited, set<size_t> &queue,
				map<size_t, vector<size_t> > &paths){
		if(visited.find(vertex) != visited.end() || vertex == size_t(-1))
			return;

//		cout << "Processing of " << vertex << endl;

		visited.insert(vertex);
		queue.erase(vertex);

		size_t rc_v = g_.IndOfRC(vertex);
		if(g_.IsVertexValid(rc_v)){
			visited.insert(rc_v);
			queue.erase(rc_v);
		}

		auto vert_for_visit = direction_definer_->GetDirectedVertices(vertex);

		if(vert_for_visit.size() == 1){
			size_t neigh = *vert_for_visit.begin();

			if(visited.find(neigh) == visited.end() && paths.find(neigh) == paths.end()){
				paths[neigh] = paths[vertex];
				paths[neigh].push_back(neigh);
				queue.insert(neigh);
			}
		}
		else
			for(auto neigh = vert_for_visit.begin(); neigh != vert_for_visit.end(); neigh++)
				if(visited.find(*neigh) == visited.end() && paths.find(*neigh) == paths.end()){
					paths[*neigh].push_back(*neigh);
					queue.insert(*neigh);
				}
	}
};

class OGD_AlternativePathProcesser : public OGD_OneVertexProcesser{
	vector<size_t> alter_path_;
	set<size_t> forbidden_vert;
	bool alter_path_is_edge;

	size_t path_start, path_end;

public:
	OGD_AlternativePathProcesser(OverlapGraph &g, OGD_DirectionDefiner *direct_definer,
			vector<size_t> alter_path) : OGD_OneVertexProcesser(g, direct_definer), alter_path_(alter_path){

		VERIFY(alter_path.size() > 1);

		for(auto e = alter_path.begin() + 1; e != alter_path.end() - 1; e++)
			forbidden_vert.insert(*e);

		alter_path_is_edge = alter_path_.size() == 2;

		path_start = *alter_path.begin(), path_end = *(alter_path.end() - 1);
	}

	void ProcessVertex(size_t vertex, set<size_t> &visited, set<size_t> &queue,
				map<size_t, vector<size_t> > &paths){
		if(visited.find(vertex) != visited.end() || vertex == size_t(-1))
			return;

		visited.insert(vertex);
		queue.erase(vertex);

		auto vert_for_visit = direction_definer_->GetDirectedVertices(vertex);

		for(auto neigh = vert_for_visit.begin(); neigh != vert_for_visit.end(); neigh++){
			if(visited.find(*neigh) == visited.end()){
				bool is_not_visit = (vertex == path_start && *neigh == path_end && alter_path_is_edge);
				if(!is_not_visit && forbidden_vert.find(*neigh) == forbidden_vert.end()){
					queue.insert(*neigh);
					paths[*neigh] = paths[vertex];
					paths[*neigh].push_back(*neigh);
				}
			}
		}
	}
};
//----------------------------------------------------------------------------------------------

class OGD_NewProcessedVertexDefiner{
protected:
	OverlapGraph &g_;

public:
	OGD_NewProcessedVertexDefiner(OverlapGraph &g) : g_(g){

	}
	virtual size_t GetNewVertex(set<size_t> &visited, set<size_t> &queue,
			map<size_t, vector<size_t> > &paths) = 0;
	virtual ~OGD_NewProcessedVertexDefiner(){

	}
};

class OGD_NewVertexInQueueDefiner : public OGD_NewProcessedVertexDefiner{

public:
	OGD_NewVertexInQueueDefiner(OverlapGraph &g) :
		OGD_NewProcessedVertexDefiner(g){
	}

	size_t GetNewVertex(set<size_t> &, set<size_t> &queue,
				map<size_t, vector<size_t> > &){
		if(queue.size() > 0)
			return *queue.begin();
		return size_t(-1);
	}
};

class OGD_SimpleNewVertexDefiner : public OGD_NewProcessedVertexDefiner{
	OGD_DirectionDefiner *direction_definer_;

public:
	OGD_SimpleNewVertexDefiner(OverlapGraph &g, OGD_DirectionDefiner *direction_definer) :
		OGD_NewProcessedVertexDefiner(g), direction_definer_(direction_definer){

	}

	size_t GetNewVertex(set<size_t> &visited, set<size_t> &queue,
				map<size_t, vector<size_t> > &paths){
		if(queue.size() > 0)
			return *queue.begin();
		else{

			auto vertices = g_.Vertices();
			for(auto v = vertices.begin(); v != vertices.end(); v++){
				if(visited.find(*v) == visited.end()){
					auto in_vertices = direction_definer_->GetAntidirectedVertices(*v);
					bool all_invisited = true;
					for(auto in_v = in_vertices.begin(); in_v != in_vertices.end(); in_v++){
						if(visited.find(*in_v) == visited.end()){
							all_invisited = false;
							break;
						}
					}

					if(all_invisited){
						paths[*v].push_back(*v);
						return *v;
					}
				}
			}

			// if vertex without antidirected edges is not exist
			// then return any unvisited vertex
			for (auto v = vertices.begin(); v != vertices.end(); v++) {
				if (visited.find(*v) == visited.end()) {
					paths[*v].push_back(*v);
					return *v;
				}
			}
		}
		return -1;
	}
};

//----------------------------------------------------------------------------------------------

class OGD_StopCondition{
protected:
	OverlapGraph &g_;
public:
	OGD_StopCondition(OverlapGraph &g) : g_(g){

	}
	virtual bool IsStop(set<size_t> &visited, set<size_t> &queue, map<size_t, vector<size_t> > &paths) = 0;
	virtual ~OGD_StopCondition(){

	}
};

class OGD_SearchedVertexIsFound : public OGD_StopCondition{
	size_t searched_vertex_;
public:
	OGD_SearchedVertexIsFound(OverlapGraph &g, size_t searched_vertex) : OGD_StopCondition(g),
		searched_vertex_(searched_vertex){

	}

	bool IsStop(set<size_t> &visited, set<size_t> &queue, map<size_t, vector<size_t> > &){
		return (visited.find(searched_vertex_) != visited.end() ||
				visited.size() == g_.VerticesCount() || queue.size() == 0);
	}
};

class OGD_NoVerticesForVisit : public OGD_StopCondition{
public:
	OGD_NoVerticesForVisit(OverlapGraph &g) : OGD_StopCondition(g){

	}

	bool IsStop(set<size_t> &visited, set<size_t> &, map<size_t, vector<size_t> > &){
		return visited.size() == g_.VerticesCount();
	}
};
//----------------------------------------------------------------------------------------------

struct OGD_Config{
	OGD_StartVerticesDefiner		* start_vert_definer;
	OGD_OneVertexProcesser 			* vertex_processer;
	OGD_NewProcessedVertexDefiner	* new_vert_definer;
	OGD_StopCondition				* stop_condition;

	OGD_Config(OGD_StartVerticesDefiner * &start_vert_definer,
			OGD_OneVertexProcesser * &vertex_processer,
			OGD_NewProcessedVertexDefiner * &new_vert_definer,
			OGD_StopCondition * &stop_condition){
		this->start_vert_definer = start_vert_definer;
		this->vertex_processer = vertex_processer;
		this->new_vert_definer = new_vert_definer;
		this->stop_condition = stop_condition;
	}

	~OGD_Config(){
		delete new_vert_definer;
		delete start_vert_definer;
		delete stop_condition;
		delete vertex_processer;
	}
};

OGD_Config CreateConfigForUniquePathsSearch(OverlapGraph &g){
	OGD_StartVerticesDefiner *start_def = new OGD_GetIsolatedAndVerticesWithoutIncoming(g);
	OGD_DirectionDefiner *direct_def = new OGD_OutgoingDirection(g);
	OGD_OneVertexProcesser *vert_proc = new OGD_UniquePathProcessing(g, direct_def);
	OGD_NewProcessedVertexDefiner *new_vert_definer = new OGD_SimpleNewVertexDefiner(g, direct_def);
	OGD_StopCondition *stop_cond = new OGD_NoVerticesForVisit(g);

	OGD_Config conf(start_def, vert_proc, new_vert_definer, stop_cond);

	return conf;
}

OGD_Config CreateContigForDijkstraFromOneVertex(OverlapGraph &g, size_t start_vertex, size_t end_vertex){
	OGD_StartVerticesDefiner *start_def = new OGD_GetParametrizedStartvertex(g, start_vertex);
	OGD_DirectionDefiner * direct_def = new OGD_OutgoingDirection(g);
	OGD_OneVertexProcesser *vert_proc = new OGD_SimpleProcessing(g, direct_def);
	OGD_NewProcessedVertexDefiner * new_vert_definer = new OGD_NewVertexInQueueDefiner(g);
	OGD_StopCondition *stop_cond = new OGD_SearchedVertexIsFound(g, end_vertex);

	OGD_Config conf(start_def, vert_proc, new_vert_definer, stop_cond);
	return conf;
}

OGD_Config CreateConfigForAlternativePathSearch(OverlapGraph &g, vector<size_t> path){

	VERIFY(path.size() > 1);
	size_t start_vertex = *(path.begin()), end_vertex = *(path.end() - 1);

	OGD_StartVerticesDefiner *start_def = new OGD_GetParametrizedStartvertex(g, start_vertex);
	OGD_DirectionDefiner * direct_def = new OGD_OutgoingDirection(g);
	OGD_OneVertexProcesser *vert_proc = new OGD_AlternativePathProcesser(g, direct_def, path);
	OGD_NewProcessedVertexDefiner * new_vert_definer = new OGD_NewVertexInQueueDefiner(g);
	OGD_StopCondition *stop_cond = new OGD_SearchedVertexIsFound(g, end_vertex);

	OGD_Config conf(start_def, vert_proc, new_vert_definer, stop_cond);
	return conf;
}

class OverlapGraphDijkstra{
	OverlapGraph &g_;
	set<size_t> visited, queue;
	map<size_t, vector<size_t> > paths;

	OGD_Config& config_;

public:
	OverlapGraphDijkstra(OverlapGraph &g, OGD_Config &config) : g_(g), config_(config){

	}

	void Run(){

//		cout << "Dijkstra run" << endl;
//		cout << "Start vertices search" << endl;
		auto start_vertices = config_.start_vert_definer->GetStartVertices();

//		cout << "Processing of start vertices" << endl;
		for(auto new_start_vertex = start_vertices.begin(); new_start_vertex != start_vertices.end();
				new_start_vertex++){
			if(visited.find(*new_start_vertex) == visited.end()){
				paths[*new_start_vertex].push_back(*new_start_vertex);
				config_.vertex_processer->ProcessVertex(*new_start_vertex, visited, queue, paths);
			}
		}

//		cout << "Dijkstra cycle starts" << endl;
		while(!config_.stop_condition->IsStop(visited, queue, paths)){
			size_t current_vertex = config_.new_vert_definer->GetNewVertex(visited, queue, paths);
			config_.vertex_processer->ProcessVertex(current_vertex, visited, queue, paths);
		}
//		cout << "Dijkstra cycle ends" << endl;
	}

	map<size_t, vector<size_t> > Paths(){
//		cout << "Paths:" << endl;
//		for(auto it = paths.begin(); it != paths.end(); it++){
//			cout << it->first << ". ";
//			auto path = it->second;
//			for(auto e = path.begin(); e != path.end(); e++)
//				cout << *e << " ";
//			cout << endl;
//		}
		return paths;
	}

	const OverlapGraph & GetGraph() { return g_; }

	~OverlapGraphDijkstra(){
	}
};

//---------------------------------------------------------------------------------

vector<vector<size_t> > DeleteRedundantEndsFromPaths(OverlapGraph &g, vector<vector<size_t> > paths){

	if(paths.size() == 0)
		return paths;

	vector<size_t> starts, ends;
	vector<bool> is_nes_start, is_nes_end;
	for(auto p = paths.begin(); p != paths.end(); p++){

		size_t cur_start = *(p->begin());
		size_t cur_end = *(p->end() - 1);

		starts.push_back(cur_start);
		ends.push_back(cur_end);

		is_nes_start.push_back(true);

		if(g.RCIdByInd(cur_start) == cur_end)
			is_nes_end.push_back(false);
		else
			is_nes_end.push_back(true);
	}

	size_t num_paths = paths.size();
	for(size_t i = 0; i < num_paths; i++){
		size_t cur_start = starts[i], cur_end = ends[i];
		for(size_t j = i + 1; j < num_paths; j++){
			size_t neig_start = starts[j], neig_end = ends[j];
			if(g.RCIdByInd(cur_start) == neig_start || g.RCIdByInd(cur_start) == neig_end)
				is_nes_start[j] = false;

			if(g.RCIdByInd(cur_end) == neig_start || g.RCIdByInd(cur_end) == neig_end)
				is_nes_end[j] = false;
		}
	}

	vector<vector<size_t> > corrected_paths;
	corrected_paths.push_back(paths[0]);

	for(size_t i = 1; i < num_paths; i++){
		if(!is_nes_start[i] || !is_nes_end[i]){
			if(paths[i].size() > 1){
				if(paths[i].size() == 2 && !is_nes_start[i] && !is_nes_end[i]){
				}
				else{
					vector<size_t> tmp;
					if(is_nes_start[i])
						tmp.push_back(paths[i][0]);
					for(size_t j = 1; j < paths[i].size() - 1; j++)
						tmp.push_back(paths[i][j]);
					if(is_nes_end[i])
						tmp.push_back(paths[i][paths[i].size() - 1]);
					corrected_paths.push_back(tmp);
				}
			}
		}
		else
			corrected_paths.push_back(paths[i]);
	}

	return corrected_paths;
}

class UniquePathsSearcher{
	OverlapGraph &g_;

	map<size_t, vector<size_t> > sh_paths;

	vector<vector<size_t> > DefineLongestPathsFromMap(){
		vector<vector<size_t> > res;
		set<size_t> used;
		while(used.size() < sh_paths.size()){
			size_t longest_path_size = 0;
			vector<size_t> longest_path;
			for(auto p = sh_paths.begin(); p != sh_paths.end(); p++){
				if(p->second.size() > longest_path_size && used.find(p->first) == used.end()){
					longest_path = p->second;
					longest_path_size = longest_path.size();
				}
			}

			for(auto v = longest_path.begin(); v != longest_path.end(); v++)
				if(sh_paths.find(*v) != sh_paths.end())
					used.insert(*v);

			res.push_back(longest_path);
		}

		return res;
	}

public:
	UniquePathsSearcher(OverlapGraph &g) : g_(g) {}

	vector<vector<size_t> > FindLongPaths(){

		OGD_Config conf = CreateConfigForUniquePathsSearch(g_);
		OverlapGraphDijkstra dijkstra(g_, conf);
		dijkstra.Run();
		sh_paths = dijkstra.Paths();

		auto long_paths = DefineLongestPathsFromMap();

		auto corrected_long_paths = DeleteRedundantEndsFromPaths(g_, long_paths);

//		cout << "Long paths" << endl;
//		for(auto p = corrected_long_paths.begin(); p != corrected_long_paths.end(); p++){
//			cout << "New path. ";
///			for(auto e = p->begin(); e != p->end(); e++)
//				cout << *e << " ";
//			cout << endl;
//		}

		return corrected_long_paths;
	}
};

class OverlapPathSearcher{
	OverlapGraph &g_;
public:
	OverlapPathSearcher(OverlapGraph &g) : g_(g) {}

	vector<size_t> GetPathAlternativeToPath(size_t start, size_t end, vector<size_t> path){
		vector<size_t> res;

		VERIFY(path.size() != 0);
		VERIFY(path[0] == start && path[path.size() - 1] == end);

		OGD_Config conf = CreateConfigForAlternativePathSearch(g_, path);
		OverlapGraphDijkstra dijkstra(g_, conf);
		dijkstra.Run();
		map<size_t, vector<size_t> > short_paths = dijkstra.Paths();

		if(short_paths.find(end) != short_paths.end()){
			res = short_paths[end];
		}

		return res;
	}

	vector<vector<size_t> > GetAlternativePaths(size_t v1, size_t v2){
		vector<vector<size_t> > paths;

//		cout << "Outgoing count - " << g_.OutgoingVerticesCount(v1) << " and incoming - " << g_.IncomingVerticesCount(v2) << endl;
		if(g_.OutgoingVerticesCount(v1) <= 1 || g_.IncomingVerticesCount(v2) <= 1)
			return paths;

		OGD_Config conf = CreateContigForDijkstraFromOneVertex(g_, v1, v2);
		OverlapGraphDijkstra dijkstra(g_, conf);
		dijkstra.Run();
		map<size_t, vector<size_t> > sh_paths = dijkstra.Paths();

		if(sh_paths.find(v2) == sh_paths.end()){
//			INFO("Path from " + ToString(v1) + " to " + ToString(v2) + " isn't found");
			return paths;
		}
		else{
			auto fst_path = sh_paths[v2];
			paths.push_back(fst_path);

			vector<size_t> snd_path = GetPathAlternativeToPath(v1, v2, fst_path);
			if(snd_path.size() != 0){
				VERIFY(snd_path[0] == v1 && snd_path[snd_path.size() - 1] == v2);
				paths.push_back(snd_path);
			}
		}

		return paths;
	}
};

//---------------------------------------------------------------------------------

void dijkstra_for_overlap_graph_test(){
	OverlapGraph g;
//	g.AddNeighVertices(1, 2, 1);
	g.AddNeighVertices(1, 3, 1);
	g.AddNeighVertices(1, 4, 1);

	g.AddNeighVertices(2, 3, 1);
	g.AddNeighVertices(4, 3, 1);

	g.AddNeighVertices(3, 4, 1);

//	OverlapPathSearcher path_searcher(g);
//	vector<int> path1;
//	path1.push_back(1);
//	path1.push_back(4);

//	auto path2 = path_searcher.GetPathAlternativeToPath(1, 4, path1);
//	for(auto v = path2.begin(); v != path2.end(); v++)
//		cout << *v << " ";
//	cout << endl;

	UniquePathsSearcher path_searcher2(g);
	auto paths = path_searcher2.FindLongPaths();

//	for(auto p = paths.begin(); p != paths.end(); p++){
//		cout << "New path. ";
//		for(auto v = p->begin(); v != p->end(); v++)
//			cout << *v << " ";
//		cout << endl;
//	}

//	auto paths_1_3 = path_searcher.GetAlternativePaths(1, 3);
//	for(auto p = paths_1_3.begin(); p != paths_1_3.end(); p++){
//		cout << "New path. ";
//		for(auto v = p->begin(); v != p->end(); v++)
//			cout << *v << " ";
//		cout << endl;
//	}
}

//---------------------------------------------------------------------------------

class OverlapGraphCorrector{
public:
	virtual size_t Correct(OverlapGraph & g) = 0;
	virtual ~OverlapGraphCorrector() { }
};


class TipClipperCorrector : OverlapGraphCorrector{
public:
	size_t Correct(OverlapGraph & g){
		auto edges = g.Edges();

		size_t deleted_edges = 0;
		for(auto e = edges.begin(); e != edges.end(); e++){
			auto start = e->first;
			auto end = e->second;

			if(g.IncomingVerticesCount(start) == 0 && g.OutgoingVerticesCount(start) == 1 &&
					g.IncomingVerticesCount(end) > 1 /*&& g.OutgoingVerticesCount(end) > 0*/){
//				cout << "Tip - " << start << " " << end << endl;
				g.DeleteVertex(start);
				deleted_edges++;
			}
			if(g.OutgoingVerticesCount(end) == 0 && g.OutgoingVerticesCount(start) > 1 &&
					/*g.IncomingVerticesCount(start) > 0 &&*/ g.IncomingVerticesCount(end) == 1){
//				cout << "Tip - " << start << " " << end << endl;
				g.DeleteVertex(end);
				deleted_edges++;
			}
		}

		return deleted_edges;
	}
};

class TransitiveReductionCorrector : public OverlapGraphCorrector{
public:
	size_t Correct(OverlapGraph & g){
		auto edges = g.Edges();
		OverlapPathSearcher ps(g);

		size_t res = 0;

		for(auto e = edges.begin(); e != edges.end(); e++){
			auto start = e->first;
			auto end = e->second;

			if(g.IsEdgeExist(start, end)){

				vector<size_t> path; path.push_back(start); path.push_back(end);
				vector<size_t> alt_path = ps.GetPathAlternativeToPath(start, end, path);

				if(alt_path.size() > 2){
					g.DeleteEdge(start, end);
					res++;
				}
			}
		}

		return res;
	}
};


class BulgeRemoverCorrector : public OverlapGraphCorrector{
public:
	size_t Correct(OverlapGraph & g){
		auto vertices = g.Vertices();
		OverlapPathSearcher ps(g);

		size_t res = 0;

		for(auto v = vertices.begin(); v != vertices.end(); v++)
			for(auto w = vertices.begin(); w != vertices.end(); w++)
				if(*v != *w && g.IsVertexExist(*v) && g.IsVertexExist(*w)){

					auto paths = ps.GetAlternativePaths(*v, *w);

					if(paths.size() > 1){

						vector<size_t> path1 = paths[0], path2 = paths[1];

						size_t w1 = 0, w2 = 0;

						for(size_t i = 0; i < path1.size() - 1; i++){
							w1 += g.GetWeightOf(path1[i], path1[i + 1]);
						}

						for(size_t i = 0; i < path2.size() - 1; i++){
							w2 += g.GetWeightOf(path2[i], path2[i + 1]);
						}

						vector<size_t> deleted_path;
						if(w1 > w2)
							deleted_path = path2;
						else
							deleted_path = path1;

						// deletion of vertices
						for(size_t i = 0; i < deleted_path.size() - 1; i++)
							g.DeleteEdge(deleted_path[i], deleted_path[i + 1]);

						// deletion of inner vertices of bulge
						for(size_t i = 1; i < deleted_path.size() - 1; i++)
							g.DeleteVertex(deleted_path[i]);

						res++;
					}
				}
		return res;
	}
};

void SimplifyOverlapGraph(OverlapGraph &overlap_graph, size_t tc_num_iter, size_t br_num_iter){

	size_t tc_res = 1, tr_res = 1;
	for(size_t i = 0; (i < tc_num_iter && (tc_res > 0 || tr_res > 0)); i++){
		TipClipperCorrector tc_corr;
		tc_res = tc_corr.Correct(overlap_graph);

		TransitiveReductionCorrector tr_corr;
		tr_res = tr_corr.Correct(overlap_graph);

		INFO(ToString(tc_res) + " tips and " + ToString(tr_res) + " transitive edges were deleted in overlap graph");
	}

	INFO("Bulge remover starts");
	BulgeRemoverCorrector br_corr;
	size_t num_bulges  = 1;
	for(size_t i = 0; (i < br_num_iter && num_bulges > 0); i++){
		num_bulges = br_corr.Correct(overlap_graph);
		INFO(ToString(num_bulges) + " bulges were deleted in overlap graph");
	}
}

}
