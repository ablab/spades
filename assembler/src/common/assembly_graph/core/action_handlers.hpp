//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __OMNI_ACTION_HANDLERS_HPP__
#define __OMNI_ACTION_HANDLERS_HPP__

#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"

#include <boost/noncopyable.hpp>
#include <string>
#include <vector>

namespace omnigraph {

/**
* ActionHandler is base listening class for graph events. All structures and information storages
* which are meant to synchronize with graph should use this structure. In order to make handler listen
* to graph events one should add it to graph listeners.
* Normally structure itself extends ActionHandler and overrides several handling methods. In
* constructor it adds itself to graph handler list and removes itself form this list in destructor.
* All events are divided into two levels: low level events and high level events.
* Low level events are addition/deletion of vertices/edges. These events should be triggered only after
* high level events when all data was already transferred and graph structure is consistent.
* High level events should be used to keep external data synchronized with graph and keep internal data
* consistent. Now high level events are merge, glue and split. This list can be extended in near future.
*/
template<typename VertexId, typename EdgeId>
class ActionHandler : private boost::noncopyable {
    const std::string handler_name_;
private:
    bool attached_;
public:
    /**
     * Create action handler with given name. With this name one can find out what tipe of handler is it.
     */
    ActionHandler(const std::string &name)
            : handler_name_(name), attached_(true) {
    }

    virtual ~ActionHandler() {
        TRACE("~ActionHandler " << handler_name_);
    }

    /**
     * Method returns name of this handler
     */
    const std::string &name() const {
        return handler_name_;
    }

    /**
     * Low level event which is triggered when vertex is added to graph.
     * @param v new vertex
     */
    virtual void HandleAdd(VertexId /*v*/) { }

    /**
     * Low level event which is triggered when edge is added to graph.
     * @param e new edge
     */
    virtual void HandleAdd(EdgeId /*e*/) { }

    /**
     * Low level event which is triggered when vertex is deleted from graph.
     * @param v vertex to delete
     */
    virtual void HandleDelete(VertexId /*v*/) { }

    /**
     * Low level event which is triggered when edge is deleted from graph.
     * @param e edge to delete
     */
    virtual void HandleDelete(EdgeId /*e*/) { }

    /**
     * High level event which is triggered when merge operation is performed on graph, which is when
     * path of edges with all inner vertices having exactly one incoming and one outgoing edge is
     * replaced with a single edge. Since this is high level operation event of creation of new edge
     * and events of deletion of old edges should not have been triggered yet when this event was triggered.
     * @param old_edges path of edges to be replaced with single edge
     * @param new_edge new edge that was added to be a replacement of path
     */
    virtual void HandleMerge(const std::vector<EdgeId> & /*old_edges*/, EdgeId /*new_edge*/) { }

    /**
     * High level event which is triggered when glue operation is performed on graph, which is when
     * edge is completely replaced with other edge. This operation is widely used in bulge removal
     * when alternative path is glued to main path. Since this is high level operation event of deletion
     * of old edge should not have been triggered yet when this event was triggered.
     * @param new_edge edge glue result
     * @param edge1 edge to be glued to edge2
     * @param edge2 edge edge1 should be glued with
     */
    virtual void HandleGlue(EdgeId /*new_edge*/, EdgeId /*edge1*/, EdgeId /*edge2*/) { }

    /**
     * High level event which is triggered when split operation is performed on graph, which is when
     * edge is split into several shorter edges. Split operation is reverse to merge operation.
     * Since this is high level operation event of deletion of old edge and events of creation of new edges
     * should not have been triggered yet when this event was triggered.
     * @param old_edge edge to be split
     * @param new_edges edges which are results of split
     */
    virtual void HandleSplit(EdgeId /*old_edge*/, EdgeId /*new_edge_1*/,
                             EdgeId /*new_edge_2*/) { }

    /**
     * Every thread safe descendant should override this method for correct concurrent graph processing.
     */
    virtual bool IsThreadSafe() const {
        return false;
    }

    bool IsAttached() const {
        return attached_;
    }

    void Attach() {
        VERIFY(!attached_);
        attached_ = true;
    }

    void Detach() {
        VERIFY(attached_);
        attached_ = false;
    }
};

template<class Graph>
class GraphActionHandler : public ActionHandler<typename Graph::VertexId,
        typename Graph::EdgeId> {
    typedef ActionHandler<typename Graph::VertexId, typename Graph::EdgeId> base;

    const Graph &g_;

protected:
    const Graph &g() const {
        return g_;
    }

public:
    GraphActionHandler(const Graph &g, const std::string &name)
            : base(name),
              g_(g) {
        TRACE("Adding new action handler: " << this->name());
        g_.AddActionHandler(this);
    }

    GraphActionHandler(const GraphActionHandler<Graph> &other)
            : base(other.name()),
              g_(other.g_) {
        TRACE("Adding new action handler: " << this->name());
        g_.AddActionHandler(this);
    }

    virtual ~GraphActionHandler() {
        TRACE("Removing action handler: " << this->name());
        if (this->IsAttached())
            this->Detach();
        g_.RemoveActionHandler(this);
    }
};

/**
* In order to support various types of graphs and make handler structure more flexible HandlerApplier
* structure was introduced. If certain implementation of graph requires special handler triggering scheme
* one can store certain extension of HandlerApplier in graph and trigger HandlerApplier methods instead
* of GraphHandler methods.
* HandlerApplier contains one method for each of graph events which define the exact way this event
* should be triggered.
*/
template<typename VertexId, typename EdgeId>
class HandlerApplier {
    typedef ActionHandler<VertexId, EdgeId> Handler;
public:

    virtual void
            ApplyAdd(Handler &handler, VertexId v) const = 0;

    virtual void
            ApplyAdd(Handler &handler, EdgeId e) const = 0;

    virtual void
            ApplyDelete(Handler &handler, VertexId v) const = 0;

    virtual void
            ApplyDelete(Handler &handler, EdgeId e) const = 0;

    virtual void ApplyMerge(Handler &handler, const std::vector<EdgeId> &old_edges,
                            EdgeId new_edge) const = 0;

    virtual void ApplyGlue(Handler &handler, EdgeId new_edge, EdgeId edge1,
                           EdgeId edge2) const = 0;

    virtual void ApplySplit(Handler &handler, EdgeId old_edge,
                            EdgeId new_edge_1, EdgeId new_edge2) const = 0;

    virtual ~HandlerApplier() {
    }
};

/**
* SimpleHandlerApplier is simple implementation of handler applier with no special filtering.
*/
template<class Graph>
class SimpleHandlerApplier : public HandlerApplier<typename Graph::VertexId,
        typename Graph::EdgeId> {
public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef ActionHandler<VertexId, EdgeId> Handler;

    void ApplyAdd(Handler &handler, VertexId v) const override {
        handler.HandleAdd(v);
    }

    void ApplyAdd(Handler &handler, EdgeId e) const override {
        handler.HandleAdd(e);
    }

    void ApplyDelete(Handler &handler, VertexId v) const override {
        handler.HandleDelete(v);
    }

    void ApplyDelete(Handler &handler, EdgeId e) const override {
        handler.HandleDelete(e);
    }

    void ApplyMerge(Handler &handler, const std::vector<EdgeId> &old_edges,
                            EdgeId new_edge) const override {
        handler.HandleMerge(old_edges, new_edge);
    }

    void ApplyGlue(Handler &handler, EdgeId new_edge, EdgeId edge1,
                           EdgeId edge2) const override {
        handler.HandleGlue(new_edge, edge1, edge2);
    }

    void ApplySplit(Handler &handler, EdgeId old_edge, EdgeId new_edge1,
                            EdgeId new_edge2) const override {
        handler.HandleSplit(old_edge, new_edge1, new_edge2);
    }

};

/**
* PairedHandlerApplier is implementation of HandlerApplier for graph with synchronization of actions
* performed with vertices/edges and its reverse-complement analogues. Thus while corresponding
* method was called only once event should be triggered twice: for the parameters with which method
* was called and for reverse-complement parameters. Also certain assertions were added for bad cases.
*/
template<class Graph>
class PairedHandlerApplier : public HandlerApplier<typename Graph::VertexId,
        typename Graph::EdgeId> {
public:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef ActionHandler<VertexId, EdgeId> Handler;

private:
    Graph &graph_;

    std::vector<EdgeId> RCPath(const std::vector<EdgeId> &path) const {
        std::vector<EdgeId> rc_path;
        rc_path.reserve(path.size());
        for (auto it = path.rbegin(), end = path.rend(); it != end; ++it) {
            rc_path.push_back(graph_.conjugate(*it));
        }
        return rc_path;
    }

public:
    PairedHandlerApplier(Graph &graph)
            : graph_(graph) {
    }

    void ApplyAdd(Handler &handler, VertexId v) const override {
        VertexId rcv = graph_.conjugate(v);
        handler.HandleAdd(v);
        if (v != rcv) {
            handler.HandleAdd(rcv);
        }
    }

    void ApplyAdd(Handler &handler, EdgeId e) const override {
        EdgeId rce = graph_.conjugate(e);
        handler.HandleAdd(e);
        if (e != rce) {
            handler.HandleAdd(rce);
        }
    }

    void ApplyDelete(Handler &handler, VertexId v) const override {
        VertexId rcv = graph_.conjugate(v);
        handler.HandleDelete(v);
        if (v != rcv) {
            handler.HandleDelete(rcv);
        }
    }

    void ApplyDelete(Handler &handler, EdgeId e) const override {
        EdgeId rce = graph_.conjugate(e);
        handler.HandleDelete(e);
        if (e != rce) {
            handler.HandleDelete(rce);
        }
    }

    void ApplyMerge(Handler &handler, const std::vector<EdgeId> &old_edges,
                            EdgeId new_edge) const override {
        EdgeId rce = graph_.conjugate(new_edge);
        handler.HandleMerge(old_edges, new_edge);
        if (new_edge != rce) {
            handler.HandleMerge(RCPath(old_edges), rce);
        }
    }

    void ApplyGlue(Handler &handler, EdgeId new_edge, EdgeId edge1,
                           EdgeId edge2) const override {
        EdgeId rc_edge1 = graph_.conjugate(edge1);
        EdgeId rc_edge2 = graph_.conjugate(edge2);
        VERIFY(edge1 != edge2);
        VERIFY(edge2 != rc_edge2);
        handler.HandleGlue(new_edge, edge1, edge2);
        if (edge1 != rc_edge1) {
            handler.HandleGlue(graph_.conjugate(new_edge), rc_edge1, rc_edge2);
        }
    }

    void ApplySplit(Handler &handler, EdgeId old_edge,
                            EdgeId new_edge_1, EdgeId new_edge2) const override {
        EdgeId rce = graph_.conjugate(old_edge);
        //VERIFY(old_edge != rce);
        handler.HandleSplit(old_edge, new_edge_1, new_edge2);
        if (old_edge != rce) {
            handler.HandleSplit(rce, graph_.conjugate(new_edge2),
                                graph_.conjugate(new_edge_1));
        }
    }

private:
    DECL_LOGGER("PairedHandlerApplier")
};

};

#endif
