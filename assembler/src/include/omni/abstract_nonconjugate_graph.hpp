//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef ABSTRACT_NONCONJUGATE_GRAPH_HPP_
#define ABSTRACT_NONCONJUGATE_GRAPH_HPP_

#include <vector>
#include <set>
#include <cstring>
#include "sequence/seq.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include "io/paired_read.hpp"
#include "omni_utils.hpp"
#include "abstract_graph.hpp"

namespace omnigraph {

template<class DataMaster>
class AbstractNonconjugateGraph;

template<class DataMaster>
class SingleEdge;

template<class DataMaster>
class SingleVertex {
private:
	typedef restricted::pure_pointer<SingleVertex<DataMaster>> VertexId;
	typedef restricted::pure_pointer<SingleEdge<DataMaster>> EdgeId;
	typedef typename DataMaster::VertexData VertexData;
public:
  typedef typename std::vector<EdgeId>::const_iterator edge_const_iterator;
private:
	friend class AbstractGraph<restricted::pure_pointer<SingleVertex<DataMaster>>, restricted::pure_pointer<SingleEdge<DataMaster>>, DataMaster> ;
	friend class AbstractNonconjugateGraph<DataMaster> ;

	vector<EdgeId> outgoing_edges_;

	vector<EdgeId> incoming_edges_;

	VertexData data_;

	size_t OutgoingEdgeCount() const {
		return outgoing_edges_.size();
	}

    edge_const_iterator out_begin() const {
        return outgoing_edges_.cbegin();
    }

    edge_const_iterator out_end() const {
        return outgoing_edges_.cend();
    }

    size_t IncomingEdgeCount() const {
        return incoming_edges_.size();
    }

    edge_const_iterator in_begin() const {
        return incoming_edges_.cbegin();
    }

    edge_const_iterator in_end() const {
        return incoming_edges_.cend();
    }

	SingleVertex(VertexData data) :
			data_(data) {
	}

	VertexData &data() {
		return data_;
	}

	void set_data(VertexData data) {
		data_ = data;
	}

	//	bool IsDeadend() {
	//		return outgoing_edges_.size() == 0;
	//	}

	void AddOutgoingEdge(EdgeId e) {
		outgoing_edges_.push_back(e);
	}

	bool RemoveOutgoingEdge(const EdgeId e) {
		auto it = outgoing_edges_.begin();
		while (it != outgoing_edges_.end() && *it != e) {
			++it;
		}
		if (it == outgoing_edges_.end()) {
			return false;
		}
		outgoing_edges_.erase(it);
		return true;
	}

	//	bool IsDeadstart() {
	//		return incoming_edges_.size() == 0;
	//	}

	void AddIncomingEdge(EdgeId e) {
		incoming_edges_.push_back(e);
	}

	bool RemoveIncomingEdge(const EdgeId e) {
		auto it = incoming_edges_.begin();
		while (it != incoming_edges_.end() && *it != e) {
			++it;
		}
		if (it == incoming_edges_.end()) {
			return false;
		}
		incoming_edges_.erase(it);
		return true;
	}
	//todo remove if nobody uses
	const vector<EdgeId> OutgoingEdgesTo(VertexId v) const {
		vector<EdgeId> result;
		for (auto it = outgoing_edges_.begin(); it != outgoing_edges_.end();
				++it) {
			if ((*it)->end() == v) {
				result.push_back(*it);
			}
		}
		return result;
	}

	~SingleVertex() {
		VERIFY(outgoing_edges_.size() == 0);
	}
};

template<class DataMaster>
class SingleEdge {
private:
	typedef restricted::pure_pointer<SingleVertex<DataMaster>> VertexId;
	typedef restricted::pure_pointer<SingleEdge<DataMaster>> EdgeId;
	typedef typename DataMaster::EdgeData EdgeData;

	friend class AbstractGraph<restricted::pure_pointer<SingleVertex<DataMaster>>, restricted::pure_pointer<SingleEdge<DataMaster>>, DataMaster> ;
	friend class AbstractNonconjugateGraph<DataMaster> ;
	//todo unfriend
	friend class SingleVertex<DataMaster> ;

	VertexId start_;
	VertexId end_;

	EdgeData data_;

	SingleEdge(VertexId start, VertexId end, const EdgeData &data) :
			start_(start), end_(end), data_(data) {
	}

	EdgeData &data() {
		return data_;
	}

	void set_data(EdgeData &data) {
		data_ = data;
	}

	void SetStartVertex(VertexId start) {
		start_ = start;
	}

	VertexId start() const {
		return start_;
	}

	void SetEndVertex(VertexId end) {
		end_ = end;
	}

	VertexId end() const {
		return end_;
	}

	~SingleEdge() {
	}
};

template<class DataMaster>
class AbstractNonconjugateGraph: public AbstractGraph<restricted::pure_pointer<SingleVertex<DataMaster>>, restricted::pure_pointer<SingleEdge<DataMaster>>, DataMaster> {
private:
	typedef AbstractGraph<restricted::pure_pointer<SingleVertex<DataMaster>>, restricted::pure_pointer<SingleEdge<DataMaster>>, DataMaster> base;
public:
	typedef typename base::VertexId VertexId;
	typedef typename base::EdgeId EdgeId;
	typedef typename base::VertexData VertexData;
	typedef typename base::EdgeData EdgeData;
	typedef typename base::VertexIterator VertexIterator;

protected:
    using base::CreateVertex;
    using base::HiddenAddEdge;

	virtual VertexId CreateVertex(const VertexData &data) {
		return VertexId(new SingleVertex<DataMaster>(data));
	}

	/*virtual */void DestroyVertex(VertexId vertex) {
		delete vertex.get();
	}

	virtual void AddVertexToGraph(VertexId vertex) {
		this->vertices_.insert(vertex);
	}

	virtual void DeleteVertexFromGraph(VertexId vertex) {
		this->vertices_.erase(vertex);
	}

	virtual VertexId HiddenAddVertex(const VertexData &data) {
		VertexId vertex = CreateVertex(data);
		AddVertexToGraph(vertex);
		return vertex;
	}

	virtual void HiddenDeleteVertex(VertexId vertex) {
		DeleteVertexFromGraph(vertex);
		DestroyVertex(vertex);
	}

    virtual void LinkIncomingEdge(VertexId v, EdgeId e) {
    	VERIFY(this->EdgeEnd(e) == VertexId(0));
    	e->SetEndVertex(v);
    	v->AddIncomingEdge(v);
    }

    virtual void LinkOutgoingEdge(VertexId v, EdgeId e) {
    	VERIFY(this->EdgeStart(e) == VertexId(0));
    	e->SetStartVertex(v);
    	v->AddOutgpoingEdge(e);
    }

	virtual EdgeId HiddenAddEdge(const EdgeData &data,
			restricted::IdDistributor * idDistributor) {
		EdgeId newEdge(new SingleEdge<DataMaster>(VertexId(0), VertexId(0), data), idDistributor);
		return newEdge;
	}

	virtual EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData &data,
			restricted::IdDistributor * idDistributor) {
		VERIFY(
				this->vertices_.find(v1) != this->vertices_.end() && this->vertices_.find(v2) != this->vertices_.end());
		EdgeId newEdge(new SingleEdge<DataMaster>(v1, v2, data), idDistributor);
		v1->AddOutgoingEdge(newEdge);
		v2->AddIncomingEdge(newEdge);
		return newEdge;
	}

	virtual void HiddenDeleteEdge(EdgeId edge) {
		VertexId start = edge->start();
		VertexId end = edge->end();
		start->RemoveOutgoingEdge(edge);
		end->RemoveIncomingEdge(edge);
		delete edge.get();
	}

	virtual vector<EdgeId> CorrectMergePath(const vector<EdgeId>& path) const {
		return path;
	}

	virtual vector<EdgeId> EdgesToDelete(const vector<EdgeId> &path) const {
		return path;
	}

	virtual vector<VertexId> VerticesToDelete(const vector<EdgeId> &path) const {
		vector<VertexId> answer;
		for (size_t i = 0; i + 1 < path.size(); i++) {
			EdgeId e = path[i + 1];
			VertexId v = EdgeStart(e);
			answer.push_back(v);
		}
		return answer;
	}


public:

	AbstractNonconjugateGraph(const DataMaster& master) :
			base(new SimpleHandlerApplier<AbstractNonconjugateGraph>(), master) {
	}

	virtual ~AbstractNonconjugateGraph() {
		//		while (!this->vertices_.empty()) {
		//			ForceDeleteVertex(*this->vertices_.begin());
		//		}
		TRACE("~AbstractNonconjugateGraph")
		for (auto it = this->SmartVertexBegin(); !it.IsEnd(); ++it) {
			ForceDeleteVertex(*it);
		}
		TRACE("~AbstractNonconjugateGraph ok")
	}

	/*virtual*/
	bool RelatedVertices(VertexId v1, VertexId v2) const {
		return v1 == v2;
	}

private:
	DECL_LOGGER("AbstractNonconjugateGraph")
};

}
#endif /* ABSTRACT_NONCONJUGATE_GRAPH_HPP_ */
