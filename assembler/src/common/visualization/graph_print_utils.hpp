//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef GRAPH_PRINTER_HPP_
#define GRAPH_PRINTER_HPP_

#include "utils/standard_base.hpp"

namespace visualization {

template<class VertexId>
struct BaseVertex {
    VertexId id_;
    string label_;
    string href_;
    string fill_color_;
    BaseVertex(VertexId id, string label, string reference, string fill_color) :id_(id), label_(label), href_(reference), fill_color_(fill_color) {
    }
};

template<class VertexId>
struct BaseEdge {
    VertexId from;
    VertexId to;
    string label;
    string color;
    BaseEdge(VertexId _from, VertexId _to, string _label, string _color) {
        from = _from;
        to = _to;
        label = _label;
        color = _color;
    }
};

class StreamRecorder {
private:
    ostream &os_;
protected:
    virtual ostream &os() {
        return os_;
    }
public:
    StreamRecorder(ostream &os) : os_(os) {
    }

    virtual ~StreamRecorder() {
    }
};

template<class Vertex, class Edge>
class GraphRecorder {
public:
    virtual void recordVertex(Vertex vertex) = 0;

    virtual void recordEdge(Edge edge) = 0;

    virtual inline void startGraphRecord(const string &name) = 0;

    virtual inline void endGraphRecord() = 0;

    virtual ~GraphRecorder(){
    }
};

template<class VertexId>
class SingleGraphRecorder : public GraphRecorder<BaseVertex<VertexId>, BaseEdge<VertexId>> {
protected:
    typedef BaseVertex<VertexId> Vertex;
    typedef BaseEdge<VertexId> Edge;
};

template<class VertexId>
class PairedGraphRecorder : public GraphRecorder<pair<BaseVertex<VertexId>, BaseVertex<VertexId>>, BaseEdge<pair<VertexId, VertexId>>> {
protected:
    typedef pair<BaseVertex<VertexId>, BaseVertex<VertexId>> Vertex;
    typedef BaseEdge<pair<VertexId, VertexId>> Edge;
};

template<class VertexId>
class DotGraphRecorder : public StreamRecorder {
public:
    DotGraphRecorder(ostream &os) : StreamRecorder(os) {
    }

protected:
    template<class vid>
    void recordVertexId(vid id) {
        this->os() << "vertex_" << id;
    }

    string IdToStr(VertexId u) {
        stringstream ss;
        ss << u;
        return ss.str();
    }

    string constructNodeId(VertexId v) {
        return constructNodePairId(v, v);
    }

    inline void recordParameter(ostream &os, const string &name, const string &value) {
        os << name << "=" << "<" << value << "> ";
    }

    inline void recordParameter(const string &name, const string &value) {
        recordParameter(this->os(), name, value);
    }

    inline void recordParameterInQuotes(ostream &os, const string &name, const string &value) {
        os << name << "=" << "\"" << value << "\" ";
    }

    inline void recordParameterInQuotes(const string &name, const string &value) {
        recordParameterInQuotes(this->os(), name, value);
    }

    inline double getColorParameter(int l, int r, double perc) {
        return l * perc + r * (1 - perc);
    }

    inline string getColor(int currentLength, int approximateLength) {
        currentLength %= approximateLength;
        int points[8][3] = {{0, 0, 1}, {0, 1, 1}, {1, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 1}, {0, 0, 1}};
        stringstream ss;
        int bound = approximateLength / 6;
        int num = currentLength / bound;
        double perc = (currentLength % bound) * 1. / bound;
        for (int i = 0; i < 3; i++) {
            ss << getColorParameter(points[num][i], points[num + 1][i], perc);
            if (i != 2)
                ss << ",";
        }
        return ss.str();
    }

};


template<class SingleVertexId>
class DotSingleGraphRecorder: public SingleGraphRecorder<SingleVertexId>, public DotGraphRecorder<SingleVertexId> {
private:
    typedef BaseVertex<SingleVertexId> Vertex;
    typedef BaseEdge<SingleVertexId> Edge;

public:
    DotSingleGraphRecorder(ostream &os) : DotGraphRecorder<SingleVertexId>(os) {
    }

    void recordVertex(Vertex vertex) {
        this->recordVertexId(vertex.id_);
        this->os() << "[";
        this->recordParameterInQuotes("label", vertex.label_);
        this->os() << ",";
        this->recordParameter("style", "filled");
        this->os() << ",";
        this->recordParameter("color", "black");
        this->os() << ",";
        if(vertex.href_ != "") {
            this->recordParameterInQuotes("href", vertex.href_);
            this->os() << ",";
        }
        this->recordParameter("fillcolor", vertex.fill_color_);
        this->os() << "]" << endl;
    }

    void recordEdge(Edge edge) {
        this->recordVertexId(edge.from);
        this->os() << "->";
        this->recordVertexId(edge.to);
        this->os() << "[";
        this->recordParameterInQuotes("label", edge.label);
        this->os() << ",";
        this->recordParameter("color", edge.color);
        this->os() << "]" << endl;
    }

    inline void startGraphRecord(const string &name) {
        this->os() << "digraph " << name << " {" << endl;
        this->os() << "node" << "[";
        this->recordParameter("fontname", "Courier");
        this->recordParameter("penwidth", "1.8");
        this->os() << "]" << endl;
    }

    inline void endGraphRecord() {
        this->os() << "}" << endl;
    }
};

template<class SingleVertexId>
class DotPairedGraphRecorder: public PairedGraphRecorder<SingleVertexId>, public DotGraphRecorder<SingleVertexId> {
private:
    typedef BaseVertex<SingleVertexId> SingleVertex;
    typedef BaseEdge<SingleVertexId> SingleEdge;
    typedef typename PairedGraphRecorder<SingleVertexId>::Vertex Vertex;
    typedef typename PairedGraphRecorder<SingleVertexId>::Edge Edge;


    string constructNodePairId(SingleVertexId u, SingleVertexId v) {
        stringstream ss;
        string u_str = this->IdToStr(u);
        string v_str = this->IdToStr(v);
        if (u == v)
            ss << u;
        else if (u_str > v_str)
            ss << v_str << "_" << u_str;
        else
            ss << u_str << "_" << v_str;
        return ss.str();
    }

    inline string constructPortCell(const string &port, string href, const string &color) {
        stringstream ss;
        ss << "<TD BORDER=\"0\" PORT = \"port_" << port << "\" ";
        this->recordParameterInQuotes(ss, "color", color);
        this->recordParameterInQuotes(ss, "bgcolor", color);
        if(href != "") {
            ss <<"href=\"" << href << "\"";
        }
        ss << "></TD>";
        return ss.str();
    }

    inline string constructLabelCell(const string &label, const string &href, const string &color) {
        stringstream ss;
        ss << "<TD BORDER=\"0\" ";
        this->recordParameterInQuotes(ss, "color", color);
        this->recordParameterInQuotes(ss, "bgcolor", color);
        if(href != "") {
            ss <<"href=\"" << href << "\"";
        }
        ss << ">"
                << label << "</TD>";
        return ss.str();
    }

    string constructComplexNodeId(string pairId, SingleVertexId v) {
        stringstream ss;
        ss << pairId << ":port_" << v;
        return ss.str();
    }

    string constructTableEntry(SingleVertex v/*, const string &label, const string &href*/) {
        stringstream ss;
        ss << "<TR>";
        ss << constructPortCell(std::to_string(v.id_) + "_in", v.href_, v.fill_color_);
        ss << constructLabelCell(v.label_, v.href_, v.fill_color_);
        ss << constructPortCell(std::to_string(v.id_) + "_out", v.href_, v.fill_color_);
        ss << "</TR>\n";
        return ss.str();
    }

    string constructReverceTableEntry(SingleVertex v/*, const string &label, const string &href*/) {
        stringstream ss;
        ss << "<TR>";
        ss << constructPortCell(std::to_string(v.id_) + "_out", v.href_, v.fill_color_);
        ss << constructLabelCell(v.label_, v.href_, v.fill_color_);
        ss << constructPortCell(std::to_string(v.id_) + "_in", v.href_, v.fill_color_);
        ss << "</TR>\n";
        return ss.str();
    }

    string constructComplexNodeLabel(Vertex v) {
        return "<TABLE BORDER=\"1\" CELLSPACING=\"0\" >\n" + constructTableEntry(v.first)
                + constructReverceTableEntry(v.second) + "</TABLE>";
    }

    string constructVertexInPairId(SingleVertexId v, SingleVertexId rc) {
        return constructComplexNodeId(constructNodePairId(v, rc), v);
    }


public:
    DotPairedGraphRecorder(ostream &os) : DotGraphRecorder<SingleVertexId>(os) {
    }

    void recordPairedVertexId(SingleVertexId id1, SingleVertexId id2) {
        this->os() << "vertex_" << constructNodePairId(id1, id2);
    }

    void recordVertex(Vertex vertex) {
        string pairLabel = constructComplexNodeLabel(vertex);
        recordPairedVertexId(vertex.first.id_, vertex.second.id_);
        this->os() << "[";
        this->recordParameter("label", constructComplexNodeLabel(vertex));
        this->os() << ",";
        this->recordParameter("color", "black");
        this->os() << ",";
        this->recordParameter("URL", "/vertex/" + std::to_string(vertex.first.id_) + ".svg");
        this->os() << "]" << endl;
    }

    void recordEdge(Edge edge) {
        this->recordVertexId(constructVertexInPairId(edge.from.first, edge.from.second));
        this->os() << "_out";
        this->os() << "->";
        this->recordVertexId(constructVertexInPairId(edge.to.first, edge.to.second));
        this->os() << "_in";
        this->os() << "[";
        this->recordParameterInQuotes("label", edge.label);
        this->os() << ",";
        this->recordParameter("color", edge.color);
        this->os() << "]" << endl;
    }

    inline void startGraphRecord(const string &name) {
        this->os() << "digraph " << name << " {" << endl;
        this->os() << "node" << "[";
        this->recordParameter("fontname", "Courier");
        this->os() << ",";
        this->recordParameter("penwidth", "1.8");
        this->os() << ",";
        this->recordParameter("shape", "plaintext");
        this->os() << "]" << endl;
    }

    inline void endGraphRecord() {
        this->os() << "}" << endl;
    }
};

}
#endif //GRAPH_PRINTER_HPP_//
