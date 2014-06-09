//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef OBSERVABLE_GRAPH_HPP_
#define OBSERVABLE_GRAPH_HPP_

#include "omni_utils.hpp"
#include "id_track_handler.hpp"
//#include <limits>

namespace omnigraph {

template<typename VertexIdT, typename EdgeIdT, typename VertexIterator/* = typename set<VertexIdT>::iterator*/>
class ObservableGraph : private boost::noncopyable {
 public:
    typedef VertexIdT VertexId;
    typedef EdgeIdT EdgeId;
    typedef VertexIterator VertexIt;
    typedef HandlerApplier<VertexId, EdgeId> Applier;
    typedef typename VertexId::type::edge_const_iterator edge_const_iterator;
    typedef SmartVertexIterator<ObservableGraph> SmartVertexIt;
    typedef SmartEdgeIterator<ObservableGraph> SmartEdgeIt;
    typedef ConstEdgeIterator<ObservableGraph> ConstEdgeIt;

    virtual void print_handlers() const {
        FOREACH (Handler* handler_ptr, action_handler_list_) {
            cout << handler_ptr->name() << endl;
        }
    }
 private:
    typedef ActionHandler<VertexId, EdgeId> Handler;

    const HandlerApplier<VertexId, EdgeId> *applier_;

    mutable vector<Handler*> action_handler_list_;

 protected:
    bool VerifyAllDetached() {
        FOREACH (Handler* handler_ptr, action_handler_list_) {
            if(handler_ptr->IsAttached()) {
            	return false;
            }
        }
        return true;
    }

 public:
 //todo make Fire* protected once again with helper friend class
    virtual void FireAddVertex(VertexId v) const {
        TRACE("FireAddVertex event of vertex inner_id=" << v.int_id() << " for " << action_handler_list_.size() << " handlers");
        FOREACH (Handler* handler_ptr, action_handler_list_) {
            if(handler_ptr->IsAttached()) {
            	TRACE("FireAddVertex to handler " << handler_ptr->name());
            	applier_->ApplyAdd(*handler_ptr, v);
            }
        }
    }

    virtual void FireAddEdge(EdgeId e) const {
        TRACE("FireAddEdge event of edge inner_id=" << e.int_id() << " for " << action_handler_list_.size() << " handlers");
        FOREACH (Handler* handler_ptr, action_handler_list_) {
            if(handler_ptr->IsAttached()) {
            	TRACE("FireAddEdge to handler " << handler_ptr->name());
            	applier_->ApplyAdd(*handler_ptr, e);
            }
        }
    }

    virtual void FireDeleteVertex(VertexId v) const {
        TRACE("FireDeleteVertex event of vertex inner_id=" << v.int_id() << " for " << action_handler_list_.size() << " handlers");
        for (auto it = action_handler_list_.rbegin(); it != action_handler_list_.rend(); ++it) {
            if((*it)->IsAttached()) {
            	TRACE("FireDeleteVertex to handler " << (*it)->name());
            	applier_->ApplyDelete(**it, v);
            }
        }
    }

    virtual void FireDeleteEdge(EdgeId e) const {
    	TRACE("FireDeleteEdge event of edge inner_id=" << e.int_id() << " for " << action_handler_list_.size() << " handlers");
        for (auto it = action_handler_list_.rbegin(); it != action_handler_list_.rend(); ++it) {
            if((*it)->IsAttached()) {
            	TRACE("FireDeleteEdge to handler " << (*it)->name());
            	applier_->ApplyDelete(**it, e);
            }
        }
        TRACE("FireDeleteEdge OK");
    }

    virtual void FireMerge(vector<EdgeId> old_edges, EdgeId new_edge) const {
        TRACE("FireMerge event, new edge inner_id=" << new_edge.int_id() << " for " << action_handler_list_.size() << " handlers");
        FOREACH (Handler* handler_ptr, action_handler_list_) {
            if(handler_ptr->IsAttached()) {
            	TRACE("FireMerge to handler " << handler_ptr->name());
            	applier_->ApplyMerge(*handler_ptr, old_edges, new_edge);
            }
        }
    }

    virtual void FireGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) const {
    	TRACE("FireGlue event, new edge inner_id=" << new_edge.int_id() << " for " << action_handler_list_.size() << " handlers");
        FOREACH (Handler* handler_ptr, action_handler_list_) {
            if(handler_ptr->IsAttached()) {
            	TRACE("FireGlue to handler " << handler_ptr->name());
            	applier_->ApplyGlue(*handler_ptr, new_edge, edge1, edge2);
            }
        }
        TRACE("FireGlue OK");
    }

    virtual void FireSplit(EdgeId edge, EdgeId new_edge1, EdgeId new_edge2) const {
        TRACE("FireSplit event, new edge1 inner_id=" << new_edge1.int_id() << ", new edge2 inner_id="
              << new_edge2.int_id() << " for " << new_edge2.int_id() << action_handler_list_.size() << " handlers");
        FOREACH (Handler* handler_ptr, action_handler_list_) {
            if(handler_ptr->IsAttached()) {
            	TRACE("FireSplit to handler " << handler_ptr->name());
            	applier_->ApplySplit(*handler_ptr, edge, new_edge1, new_edge2);
            }
        }
    }

    ObservableGraph(HandlerApplier<VertexId, EdgeId> *applier)
            : applier_(applier)/*, element_order_(*this)*/{
    }

    virtual ~ObservableGraph() {
        TRACE("~ObservableGraph")
        delete applier_;
        TRACE("~ObservableGraph ok")
    }

    void AddActionHandler(Handler* action_handler) const {
#pragma omp critical(action_handler_list_modification)
        {
            TRACE("Action handler " << action_handler->name() << " added");
            if (find(action_handler_list_.begin(), action_handler_list_.end(),
                     action_handler) != action_handler_list_.end()) {
                VERIFY_MSG(
                        false,
                        "Action handler " << action_handler->name() << " has already been added");
            } else {
                action_handler_list_.push_back(action_handler);
            }
        }
    }

    bool RemoveActionHandler(const Handler* action_handler) const {
        bool result = false;
#pragma omp critical(action_handler_list_modification)
        {
            auto it = std::find(action_handler_list_.begin(),
                                action_handler_list_.end(), action_handler);
            if (it != action_handler_list_.end()) {
                action_handler_list_.erase(it);
                TRACE("Action handler " << action_handler->name() << " removed");
                result = true;
            } else {
                TRACE("Action handler " << action_handler->name() << " wasn't found among graph action handlers");
            }
        }

        return result;
    }

    virtual VertexIterator begin() const = 0;

    virtual VertexIterator end() const = 0;

    virtual edge_const_iterator out_begin(VertexId v) const = 0;

    virtual edge_const_iterator out_end(VertexId v) const = 0;

    template<typename Comparator>
    SmartVertexIterator<ObservableGraph, Comparator> SmartVertexBegin(
            const Comparator& comparator) const {
        return SmartVertexIterator<ObservableGraph, Comparator>(*this,
                                                                comparator);
    }

    SmartVertexIterator<ObservableGraph> SmartVertexBegin() const {
        return SmartVertexIterator<ObservableGraph>(*this);
    }

    template<typename Comparator>
    SmartEdgeIterator<ObservableGraph, Comparator> SmartEdgeBegin(
            const Comparator& comparator) const {
        return SmartEdgeIterator<ObservableGraph, Comparator>(*this, comparator);
    }

    SmartEdgeIterator<ObservableGraph> SmartEdgeBegin() const {
        return SmartEdgeIterator<ObservableGraph>(*this);
    }

    ConstEdgeIterator<ObservableGraph> ConstEdgeBegin() const {
        return ConstEdgeIterator<ObservableGraph>(*this);
    }

    //Use very carefully!
    //todo remove
    void FireProject(EdgeId edge1, EdgeId edge2) {
        FireGlue(edge2, edge1, edge2);
    }

    const Applier& GetHandlerApplier() const {
        return *applier_;
    }

    bool AllHandlersThreadSafe() const {
        BOOST_FOREACH(Handler* handler, action_handler_list_) {
            if (handler->IsAttached() && !handler->IsThreadSafe()) {
                return false;
            }
        }
        return true;
    }

    // TODO: for debug. remove.
    void PrintHandlersNames() const {
        BOOST_FOREACH(Handler* handler, action_handler_list_) {
            cout << handler->name() << " attached=" << handler->IsAttached() << endl;
        }
    }

 private:
    DECL_LOGGER("ObservableGraph")
};
}
#endif /* OBSERVABLE_GRAPH_HPP_ */
