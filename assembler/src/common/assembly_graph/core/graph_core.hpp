//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <set>
#include "utils/verify.hpp"
#include "utils/logger/logger.hpp"
#include "order_and_law.hpp"
#include <boost/iterator/iterator_facade.hpp>
#include "utils/simple_tools.hpp"

namespace omnigraph {

using std::vector;
template<class DataMaster>
class GraphCore;

template<class DataMaster>
class ConstructionHelper;

template<class T>
class PairedElementManipulationHelper;

template<class DataMaster>
class PairedVertex;

template<class DataMaster>
class PairedEdge;

template<class DataMaster>
class PairedEdge {
 private:
    typedef typename DataMaster::EdgeData EdgeData;
    typedef restricted::pure_pointer<PairedEdge<DataMaster>> EdgeId;
    typedef restricted::pure_pointer<PairedVertex<DataMaster>> VertexId;
    friend class GraphCore<DataMaster>;
    friend class ConstructionHelper<DataMaster>;
    friend class PairedElementManipulationHelper<EdgeId>;
    //todo unfriend
    friend class PairedVertex<DataMaster>;
    VertexId end_;
    EdgeData data_;
    EdgeId conjugate_;

    PairedEdge(VertexId end, const EdgeData &data)
            : end_(end),
              data_(data) {
    }

    EdgeData &data() {
        return data_;
    }

    void set_data(const EdgeData &data) {
        data_ = data;
    }

    VertexId end() const {
        return end_;
    }

    VertexId start() const {
        return conjugate_->end()->conjugate();
    }

    void set_conjugate(EdgeId conjugate) {
        conjugate_ = conjugate;
    }

    void SetEndVertex(VertexId end) {
        end_ = end;
    }

public:
    EdgeId conjugate() const {
        return conjugate_;
    }

    size_t length(size_t k) const {
        return data_.size() - k;
    }
};

template<class DataMaster>
class PairedVertex {
private:
    typedef typename DataMaster::VertexData VertexData;
    typedef restricted::pure_pointer<PairedEdge<DataMaster>> EdgeId;
    typedef restricted::pure_pointer<PairedVertex<DataMaster>> VertexId;
    typedef typename std::vector<EdgeId>::const_iterator edge_raw_iterator;

    class conjugate_iterator : public boost::iterator_facade<conjugate_iterator,
            EdgeId, boost::forward_traversal_tag, EdgeId> {
    public:
        explicit conjugate_iterator(edge_raw_iterator it,
                                    bool conjugate = false)
                : it_(it),
                  conjugate_(conjugate) {
        }

        //todo do we need it?
        conjugate_iterator()
                : conjugate_(false) {
        }

    private:
        friend class boost::iterator_core_access;

        void increment() {
            it_++;
        }

        bool equal(const conjugate_iterator &other) const {
            return other.it_ == it_ && other.conjugate_ == conjugate_;
        }

        EdgeId dereference() const {
            return (conjugate_ ? (*it_)->conjugate() : *it_);
        }

        edge_raw_iterator it_;
        bool conjugate_;
    };

public:
    typedef conjugate_iterator edge_const_iterator;

private:
    friend class GraphCore<DataMaster>;
    friend class ConstructionHelper<DataMaster>;
    friend class PairedEdge<DataMaster>;
    friend class PairedElementManipulationHelper<VertexId>;
    friend class conjugate_iterator;

    std::vector<EdgeId> outgoing_edges_;

    VertexId conjugate_;

    VertexData data_;

    bool IsMinimal() const {
        return conjugate_->conjugate_ <= conjugate_;
    }

    VertexId conjugate() const {
        return conjugate_;
    }

    void set_conjugate(VertexId conjugate) {
        conjugate_ = conjugate;
    }

    size_t OutgoingEdgeCount() const {
        return outgoing_edges_.size();
    }

    edge_const_iterator out_begin() const {
        return edge_const_iterator(outgoing_edges_.cbegin(), false);
    }

    edge_const_iterator out_end() const {
        return edge_const_iterator(outgoing_edges_.cend(), false);
    }

    size_t IncomingEdgeCount() const {
        return conjugate_->OutgoingEdgeCount();
    }

    size_t IncomingEdgesCount() const {
        return conjugate_->OutgoingEdgeCount();
    }

    edge_const_iterator in_begin() const {
        return edge_const_iterator(conjugate_->outgoing_edges_.cbegin(), true);
    }

    edge_const_iterator in_end() const {
        return edge_const_iterator(conjugate_->outgoing_edges_.cend(), true);
    }

    PairedVertex(VertexData data)
            : data_(data) {
    }

    VertexData &data() {
        return data_;
    }

    void set_data(VertexData data) {
        data_ = data;
    }

    const std::vector<EdgeId> OutgoingEdgesTo(VertexId v) const {
        vector<EdgeId> result;
        for (auto it = outgoing_edges_.begin(); it != outgoing_edges_.end(); ++it) {
            if ((*it)->end() == v) {
                result.push_back(*it);
            }
        }
        return result;
    }

    void AddOutgoingEdge(EdgeId e) {
        outgoing_edges_.insert(std::upper_bound(outgoing_edges_.begin(), outgoing_edges_.end(), e), e);
        //outgoing_edges_.push_back(e);
    }

    bool RemoveOutgoingEdge(const EdgeId e) {
        auto it = std::find(outgoing_edges_.begin(), outgoing_edges_.end(), e);
        if (it == outgoing_edges_.end())
            return false;

        outgoing_edges_.erase(it);
        return true;
    }

    ~PairedVertex() {
        VERIFY(outgoing_edges_.size() == 0);
    }
};

template<class DataMaster>
class GraphCore: private boost::noncopyable {
public:
    typedef DataMaster DataMasterT;
    typedef typename DataMasterT::VertexData VertexData;
    typedef typename DataMasterT::EdgeData EdgeData;
    typedef restricted::pure_pointer<PairedEdge<DataMaster>> EdgeId;
    typedef restricted::pure_pointer<PairedVertex<DataMaster>> VertexId;
    typedef typename std::set<VertexId>::const_iterator VertexIt;
    typedef typename PairedVertex<DataMaster>::edge_const_iterator edge_const_iterator;

private:
   restricted::LocalIdDistributor id_distributor_;
   DataMaster master_;
   std::set<VertexId> vertices_;

   friend class ConstructionHelper<DataMaster>;
public:
   VertexIt begin() const {
       return vertices_.begin();
   }

   VertexIt end() const {
       return vertices_.end();
   }

   const std::set<VertexId>& vertices() const {
       return vertices_;
   }

   size_t size() const {
       return vertices_.size();
   }

   edge_const_iterator out_begin(VertexId v) const {
       return v->out_begin();
   }

   edge_const_iterator out_end(VertexId v) const {
       return v->out_end();
   }

   edge_const_iterator in_begin(VertexId v) const {
       return v->in_begin();
   }

   edge_const_iterator in_end(VertexId v) const {
       return v->in_end();
   }

private:
   void DeleteVertexFromGraph(VertexId vertex) {
       this->vertices_.erase(vertex);
       this->vertices_.erase(conjugate(vertex));
   }

   void DestroyVertex(VertexId vertex) {
       VertexId conjugate = vertex->conjugate();
       delete vertex.get();
       delete conjugate.get();
   }

   bool AdditionalCompressCondition(VertexId v) const {
       return !(EdgeEnd(GetUniqueOutgoingEdge(v)) == conjugate(v) && EdgeStart(GetUniqueIncomingEdge(v)) == conjugate(v));
   }

protected:

   VertexId CreateVertex(const VertexData& data1, const VertexData& data2, restricted::IdDistributor& id_distributor) {
       VertexId vertex1(new PairedVertex<DataMaster>(data1), id_distributor);
       VertexId vertex2(new PairedVertex<DataMaster>(data2), id_distributor);
       vertex1->set_conjugate(vertex2);
       vertex2->set_conjugate(vertex1);
       return vertex1;
   }

   VertexId CreateVertex(const VertexData &data, restricted::IdDistributor &id_distributor) {
       return CreateVertex(data, master_.conjugate(data), id_distributor);
   }

    VertexId CreateVertex(const VertexData &data) {
        return CreateVertex(data, id_distributor_);
    }

    void AddVertexToGraph(VertexId vertex) {
        vertices_.insert(vertex);
        vertices_.insert(conjugate(vertex));
    }

    VertexId HiddenAddVertex(const VertexData& data, restricted::IdDistributor& id_distributor) {
        VertexId vertex = CreateVertex(data, id_distributor);
        AddVertexToGraph(vertex);
        return vertex;
    }

    VertexId HiddenAddVertex(const VertexData& data) {
        return HiddenAddVertex(data, id_distributor_);
    }

    void HiddenDeleteVertex(VertexId vertex) {
        DeleteVertexFromGraph(vertex);
        DestroyVertex(vertex);
    }

    /////////////////////////low-level ops (move to helper?!)

    ////what with this method?
    EdgeId AddSingleEdge(VertexId v1, VertexId v2, const EdgeData &data,
                         restricted::IdDistributor &idDistributor) {
        EdgeId newEdge(new PairedEdge<DataMaster>(v2, data), idDistributor);
        if (v1 != VertexId(0))
            v1->AddOutgoingEdge(newEdge);
        return newEdge;
    }

    EdgeId HiddenAddEdge(const EdgeData& data, restricted::IdDistributor& id_distributor) {
        EdgeId result = AddSingleEdge(VertexId(0), VertexId(0), data, id_distributor);
        if (this->master().isSelfConjugate(data)) {
            result->set_conjugate(result);
            return result;
        }
        EdgeId rcEdge = AddSingleEdge(VertexId(0), VertexId(0), this->master().conjugate(data), id_distributor);
        result->set_conjugate(rcEdge);
        rcEdge->set_conjugate(result);
        return result;
    }

    EdgeId HiddenAddEdge(const EdgeData &data) {
        return HiddenAddEdge(data, id_distributor_);
    }

    EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData& data, restricted::IdDistributor& id_distributor) {
        //      todo was suppressed for concurrent execution reasons (see concurrent_graph_component.hpp)
        //      VERIFY(this->vertices_.find(v1) != this->vertices_.end() && this->vertices_.find(v2) != this->vertices_.end());
        EdgeId result = AddSingleEdge(v1, v2, data, id_distributor);
        if (this->master().isSelfConjugate(data) && (v1 == conjugate(v2))) {
            //              todo why was it removed???
            //          Because of some split issues: when self-conjugate edge is split armageddon happends
            //          VERIFY(v1 == conjugate(v2));
            //          VERIFY(v1 == conjugate(v2));
            result->set_conjugate(result);
            return result;
        }
        EdgeId rcEdge = AddSingleEdge(v2->conjugate(), v1->conjugate(), this->master().conjugate(data), id_distributor);
        result->set_conjugate(rcEdge);
        rcEdge->set_conjugate(result);
        return result;
    }

    EdgeId HiddenAddEdge(VertexId v1, VertexId v2, const EdgeData &data) {
        return HiddenAddEdge(v1, v2, data, id_distributor_);
    }

    void HiddenDeleteEdge(EdgeId edge) {
        TRACE("Hidden delete edge " << edge.int_id());
        EdgeId rcEdge = conjugate(edge);
        VertexId rcStart = conjugate(edge->end());
        VertexId start = conjugate(rcEdge->end());
        start->RemoveOutgoingEdge(edge);
        rcStart->RemoveOutgoingEdge(rcEdge);
        if (edge != rcEdge) {
            delete rcEdge.get();
        }
        delete edge.get();
    }

    void HiddenDeletePath(const std::vector<EdgeId>& edgesToDelete, const std::vector<VertexId>& verticesToDelete) {
        for (auto it = edgesToDelete.begin(); it != edgesToDelete.end(); ++it)
            HiddenDeleteEdge(*it);
        for (auto it = verticesToDelete.begin(); it != verticesToDelete.end(); ++it)
            HiddenDeleteVertex(*it);
    }

public:

    GraphCore(const DataMaster& master) : master_(master) {
    }

    virtual ~GraphCore() {
        VERIFY(size() == 0);
    }

    class IteratorContainer {
    public:
        typedef edge_const_iterator const_iterator;
    private:
        const_iterator begin_;
        const_iterator end_;
    public:
        IteratorContainer(const_iterator begin, const_iterator end) :
            begin_(begin), end_(end) {

        }

        const_iterator begin() const {
            return begin_;
        }

        const_iterator end() const {
            return end_;
        }
    };

    restricted::LocalIdDistributor &GetGraphIdDistributor() {
        return id_distributor_;
    }

    const restricted::LocalIdDistributor &GetGraphIdDistributor() const {
        return id_distributor_;
    }

    size_t int_id(EdgeId edge) const {
        return edge.int_id();
    }

    size_t int_id(VertexId vertex) const {
        return vertex.int_id();
    }

    const DataMaster& master() const {
        return master_;
    }

    const EdgeData& data(EdgeId edge) const {
        return edge->data();
    }

    const VertexData& data(VertexId v) const {
        return v->data();
    }

    EdgeData& data(EdgeId edge) {
        return edge->data();
    }

    VertexData& data(VertexId v) {
        return v->data();
    }

    size_t OutgoingEdgeCount(VertexId v) const {
        return v->OutgoingEdgeCount();
    }

    IteratorContainer OutgoingEdges(VertexId v) const {
        //INFO("Outgoing");
        return IteratorContainer(out_begin(v), out_end(v));
    }

    size_t IncomingEdgeCount(VertexId v) const {
        return v->IncomingEdgeCount();
    }

    IteratorContainer IncomingEdges(VertexId v) const {
        return IteratorContainer(in_begin(v), in_end(v));
    }

    std::vector<EdgeId> GetEdgesBetween(VertexId v, VertexId u) const {
        return v->OutgoingEdgesTo(u);
    }

    bool RelatedVertices(VertexId v1, VertexId v2) const {
        return v1 == v2 || v1 == conjugate(v2);
    }

    ////////////////////////edge information
    VertexId EdgeStart(EdgeId edge) const {
        return edge->start();
    }

    VertexId EdgeEnd(EdgeId edge) const {
        //INFO("Edge end");
        return edge->end();
    }

    VertexId conjugate(VertexId v) const {
        return v->conjugate();
    }

    EdgeId conjugate(EdgeId edge) const {
        return edge->conjugate();
    }

    size_t length(const EdgeId edge) const {
        return master_.length(data(edge));
    }

    size_t length(const VertexId v) const {
        return master_.length(data(v));
    }

    //////////////////////shortcut methods

    std::vector<EdgeId> IncidentEdges(VertexId v) const {
        vector<EdgeId> answer;
        push_back_all(answer, IncomingEdges(v));
        push_back_all(answer, OutgoingEdges(v));
        return answer;
    }

    EdgeId GetUniqueOutgoingEdge(VertexId v) const {
        VERIFY(CheckUniqueOutgoingEdge(v));
        return *out_begin(v);
    }

    bool CheckUniqueIncomingEdge(VertexId v) const {
        return IncomingEdgeCount(v) == 1;
    }

    EdgeId GetUniqueIncomingEdge(VertexId v) const {
        VERIFY(CheckUniqueIncomingEdge(v));
        return *in_begin(v);
    }

    bool CheckUniqueOutgoingEdge(VertexId v) const {
        return OutgoingEdgeCount(v) == 1;
    }

    bool IsDeadEnd(VertexId v) const {
        return OutgoingEdgeCount(v) == 0;
    }

    bool IsDeadStart(VertexId v) const {
        return IncomingEdgeCount(v) == 0;
    }
    
    bool CanCompressVertex(VertexId v) const {
        //      TRACE("Compress vertex check: ");
        //      TRACE("Outgoing check: " << (OutgoingEdgeCount(v) == 1));
        //      TRACE("Outgoing check: " << (CheckUniqueOutgoingEdge(v)));
        //      TRACE("Incoming check: " << (IncomingEdgeCount(v) == 1));
        //      TRACE("Incoming check: " << (CheckUniqueIncomingEdge(v) == 1));
        //      if((OutgoingEdgeCount(v) == 1) && (IncomingEdgeCount(v) == 1)) {
        //          TRACE("Loop check: " << (GetUniqueOutgoingEdge(v) != GetUniqueIncomingEdge(v)));
        //          TRACE("Additional check: " << AdditionalCompressCondition(v));
        //      }
        return OutgoingEdgeCount(v) == 1 && IncomingEdgeCount(v) == 1 && 
            GetUniqueOutgoingEdge(v) != GetUniqueIncomingEdge(v) &&
            AdditionalCompressCondition(v);
    }
    
    //////////////////////printing
    std::string str(const EdgeId e) const {
//      return master_.str(data(edge));
        std::stringstream ss;
        ss << int_id(e) << " (" << length(e) << ")";
        return ss.str();
    }

    std::string str(const VertexId v) const {
//      return master_.str(data(v));
        return ToString(int_id(v));
    }

    std::string detailed_str(const VertexId v) const {
        std::stringstream ss;
        ss << str(v) << ";";
        ss << "Incoming edges" << str(IncomingEdges(v)) << "; ";
        ss << "Outgoing edges" << str(OutgoingEdges(v)) << ";";
        return ss.str();
    }

    std::string detailed_str(const std::vector<EdgeId>& path) const {
        std::stringstream ss;
        ss << "Path: ";
        ss << "Vertex " << detailed_str(EdgeStart(path[0])) << " | ";
        for (auto it = path.begin(); it != path.end(); ++it) {
            EdgeId e = *it;
            ss << "Edge " << str(e) << " | ";
            ss << "Vertex " << detailed_str(EdgeEnd(e)) << " | ";
        }
        return ss.str();
    }

    template<class Container>
    std::string str(const Container& container) const {
        return str(container.begin(), container.end());
    }

    template<class It>
    std::string str(It begin, It end) const {
        std::stringstream ss;
        std::string delim = "";
        for (auto it = begin; it != end; ++it) {
            ss << delim << str(*it);
            delim = ", ";
        }
        return ss.str();
    }

private:
    DECL_LOGGER("GraphCore");
};

}
