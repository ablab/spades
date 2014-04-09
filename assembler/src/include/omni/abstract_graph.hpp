//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "abstract_editable_graph.hpp"
#include "id_track_handler.hpp"

namespace omnigraph {

template<typename VertexIdT, typename EdgeIdT, class DataMasterT>
class AbstractGraph : public AbstractEditableGraph<VertexIdT, EdgeIdT,
        DataMasterT, typename set<VertexIdT>::const_iterator> {
    typedef AbstractEditableGraph<VertexIdT, EdgeIdT, DataMasterT,
            typename set<VertexIdT>::const_iterator> base;

 public:
    typedef VertexIdT VertexId;
    typedef EdgeIdT EdgeId;
    typedef DataMasterT DataMaster;
    typedef typename DataMaster::VertexData VertexData;
    typedef typename DataMaster::EdgeData EdgeData;
    typedef typename base::VertexIterator VertexIterator;
    typedef typename VertexId::type::edge_const_iterator edge_const_iterator;
    using base::str;

 protected:

    set<VertexId> vertices_;

 public:
    AbstractGraph(HandlerApplier<VertexId, EdgeId>* applier,
                  const DataMaster& master)
            : base(applier, master) {
    }

    virtual ~AbstractGraph() {
        TRACE("~AbstractGraph");
    }

    virtual edge_const_iterator out_begin(VertexId v) const {
        return v->out_begin();
    }

    virtual edge_const_iterator out_end(VertexId v) const {
        return v->out_end();
    }

    virtual edge_const_iterator in_begin(VertexId v) const {
        return v->in_begin();
    }

    virtual edge_const_iterator in_end(VertexId v) const {
        return v->in_end();
    }

    virtual size_t OutgoingEdgeCount(VertexId v) const {
        return v->OutgoingEdgeCount();
    }

    virtual size_t IncomingEdgeCount(VertexId v) const {
        return v->IncomingEdgeCount();
    }

    virtual vector<EdgeId> GetEdgesBetween(VertexId v, VertexId u) const {
        return v->OutgoingEdgesTo(u);
    }

    virtual const EdgeData& data(EdgeId edge) const {
        return edge->data();
    }

    virtual const VertexData& data(VertexId v) const {
        return v->data();
    }

    virtual EdgeData& data(EdgeId edge) {
        return edge->data();
    }

    virtual VertexData& data(VertexId v) {
        return v->data();
    }

    virtual VertexIterator begin() const {
        return vertices_.begin();
    }

    virtual VertexIterator end() const {
        return vertices_.end();
    }

    const set<VertexId> &vertices() const {
	    return vertices_;
	}

	size_t size() const {
		return vertices_.size();
	}

    virtual VertexId EdgeStart(EdgeId edge) const {
        return edge->start();
    }

    virtual VertexId EdgeEnd(EdgeId edge) const {
        //INFO("Edge end");
        return edge->end();
    }

 protected:

    virtual void AddVertexToGraph(VertexId vertex) = 0;

    virtual void DeleteVertexFromGraph(VertexId vertex) = 0;

 private:
    DECL_LOGGER("AbstractGraph")
    ;
};
}
