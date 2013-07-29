// Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: David Alexander
// (Based on the original "Partial Order Aligner" by Lee, Grasso, and Sharlow,
//  and an implementation in C# by Patrick Marks)

#include "Poa/PoaGraph.hpp"

#include <boost/config.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>
#include <cassert>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

#include "Poa/PoaConfig.hpp"
#include "Types.hpp"
#include "Utils.hpp"
#include "Mutation.hpp"

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::endl;

///
///  Boost graph library typedefs, properties, and graphviz output.
///
using namespace boost; // NOLINT

namespace boost
{
    enum vertex_info_t { vertex_info = 424 };  // a unique #
    BOOST_INSTALL_PROPERTY(vertex, info);
}



namespace ConsensusCore
{
    struct PoaNode
    {
        char Base;
        int Reads;
        int SpanningReads;
        float Score;
        float ReachingScore;
        bool IsInConsensus;

        void Init(char base, int reads)
        {
            this->Base = base;
            this->Reads = reads;
            this->SpanningReads = 0;
            this->Score = 0;
            this->ReachingScore = 0;
            this->IsInConsensus = false;
        }

        explicit PoaNode(char base)
        {
            Init(base, 1);
        }

        PoaNode(char base, int reads)
        {
            Init(base, reads);
        }
    };

    typedef adjacency_list < setS, vecS, bidirectionalS,
                             property<vertex_info_t, PoaNode*> > BoostGraph;
    typedef graph_traits<BoostGraph>::edge_descriptor Edge;
    typedef graph_traits<BoostGraph>::vertex_descriptor Vertex;
    typedef property_map<BoostGraph, vertex_info_t>::type VertexInfoMap;
    static const Vertex null_vertex = graph_traits<BoostGraph>::null_vertex();
}


namespace boost
{
    using ConsensusCore::VertexInfoMap;
    using boost::format;

    class my_label_writer
    {
    public:
        my_label_writer(VertexInfoMap map, bool color, bool verbose)
            : map_(map),
              color_(color),
              verbose_(verbose)
        {}

        template <class VertexOrEdge>
        void operator()(std::ostream& out, const VertexOrEdge& v) const
        {
            std::string nodeColoringAttribute =
                (color_ && map_[v]->IsInConsensus?
                 " style=\"filled\", fillcolor=\"lightblue\" ," : "");
            if (!verbose_)
            {
                out << format("[shape=Mrecord,%s label=\"{ %c | %d }\"]")
                    % nodeColoringAttribute
                    % map_[v]->Base
                    % map_[v]->Reads;
            }
            else
            {
                out <<  format("[shape=Mrecord,%s label=\"{ "
                               "{ %d | %c } |"
                               "{ %d | %d } |"
                               "{ %0.2f | %0.2f } }\"]")
                    % nodeColoringAttribute
                    % v % map_[v]->Base
                    % map_[v]->Reads % map_[v]->SpanningReads
                    % map_[v]->Score % map_[v]->ReachingScore;
            }
        }
    private:
        VertexInfoMap map_;
        bool color_;
        bool verbose_;
    };
}


namespace ConsensusCore
{
    enum MoveType
    {
        InvalidMove,  // Invalid move reaching ^ (start)
        StartMove,    // Start move: ^ -> vertex in row 0 of local alignment
        EndMove,      // End move: vertex -> $ in row 0 of local alignment, or
                      //  in global alignment, terminal vertex -> $
        MatchMove,
        MismatchMove,
        DeleteMove,
        ExtraMove
    };

    struct AlignmentColumn : noncopyable
    {
        Vertex CurrentVertex;
        vector<float> Score;
        vector<MoveType> ReachingMove;
        vector<Vertex> PreviousVertex;

        AlignmentColumn(Vertex vertex, int len)
            : CurrentVertex(vertex),
              Score(len, -FLT_MAX),
              ReachingMove(len, InvalidMove),
              PreviousVertex(len, null_vertex)
        {}

        ~AlignmentColumn()
        {}
    };

    typedef unordered_map<Vertex, const AlignmentColumn*> AlignmentColumnMap;

    //
    // Graph::Impl methods
    //

    class PoaGraph::Impl
    {
        BoostGraph g_;
        VertexInfoMap vertexInfoMap_;
        Vertex enterVertex_;
        Vertex exitVertex_;
        vector<std::string> sequences_;

        void repCheck();

        //
        // utility routines
        //
        const AlignmentColumn*
        makeAlignmentColumn(Vertex v,
                            const AlignmentColumnMap& alignmentColumnForVertex,
                            const std::string& sequence,
                            const PoaConfig& config);

        const AlignmentColumn*
        makeAlignmentColumnForExit(Vertex v,
                                   const AlignmentColumnMap& alignmentColumnForVertex,
                                   const std::string& sequence,
                                   const PoaConfig& config);

    public:
        Impl();
        ~Impl();
        void AddSequence(const std::string& sequence, const PoaConfig& config);

        // TODO(dalexander): make this const
        tuple<string, float, vector< pair<Mutation*, float> >*>
        FindConsensus(const PoaConfig& config);

        int NumSequences() const;
        string ToGraphViz(int flags) const;
        void WriteGraphVizFile(string filename, int flags) const;
    };


    PoaGraph::Impl::Impl()
    {
        vertexInfoMap_ = get(vertex_info, g_);
        enterVertex_ = add_vertex(g_);
        vertexInfoMap_[enterVertex_] = new PoaNode('^', 0);
        exitVertex_ = add_vertex(g_);
        vertexInfoMap_[exitVertex_] = new PoaNode('$', 0);
    }

    PoaGraph::Impl::~Impl()
    {
        foreach (Vertex v, vertices(g_))
        {
            delete vertexInfoMap_[v];
        }
    }

    void PoaGraph::Impl::repCheck()
    {
        // assert the representation invariant for the object
        foreach (Vertex v, vertices(g_))
        {
            if (v == enterVertex_)
            {
                assert(in_degree(v, g_) == 0);
                assert(out_degree(v, g_) > 0 || NumSequences() == 0);
            }
            else if (v == exitVertex_)
            {
                assert(in_degree(v, g_) > 0 || NumSequences() == 0);
                assert(out_degree(v, g_) == 0);
            }
            else
            {
                assert(in_degree(v, g_) > 0);
                assert(out_degree(v, g_) > 0);
            }
            assert(vertexInfoMap_[v] != NULL);
        }
    }


    static inline vector<const AlignmentColumn*>
    getPredecessorColumns(BoostGraph& g,
                          Vertex v,
                          const AlignmentColumnMap& alignmentColumnForVertex)
    {
        vector<const AlignmentColumn*> predecessorColumns;
        const AlignmentColumn* predCol;
        foreach (Edge e, in_edges(v, g))
        {
            Vertex u = source(e, g);
            predCol = alignmentColumnForVertex.at(u);
            assert(predCol != NULL);
            predecessorColumns.push_back(predCol);
        }
        return predecessorColumns;
    }

    const AlignmentColumn*
    PoaGraph::Impl::makeAlignmentColumnForExit(Vertex v,
                                               const AlignmentColumnMap& alignmentColumnForVertex,
                                               const std::string& sequence,
                                               const PoaConfig& config)
    {
        assert(out_degree(v, g_) == 0);

        // this is kind of unnecessary as we are only actually using one entry in this column
        int I = sequence.length();
        AlignmentColumn* curCol = new AlignmentColumn(v, I + 1);

        float bestScore = -FLT_MAX;
        Vertex prevVertex = null_vertex;

        // Under local alignment the vertex $ can be reached from
        // any other vertex in one step via the End move--not just its
        // predecessors in the graph
        if (config.UseLocalAlignment)
        {
            foreach (Vertex u, vertices(g_))
            {
                if (u != exitVertex_)
                {
                    const AlignmentColumn* predCol = alignmentColumnForVertex.at(u);
                    if (predCol->Score[I] > bestScore)
                    {
                        bestScore = predCol->Score[I];
                        prevVertex = predCol->CurrentVertex;
                    }
                }
            }
        }
        else
        {
            // regular predecessors
            vector<const AlignmentColumn*> predecessorColumns  =
                    getPredecessorColumns(g_, v, alignmentColumnForVertex);
            foreach (const AlignmentColumn * predCol, predecessorColumns)
            {
                if (predCol->Score[I] > bestScore)
                {
                    bestScore = predCol->Score[I];
                    prevVertex = predCol->CurrentVertex;
                }
            }
        }
        assert(prevVertex != null_vertex);
        curCol->Score[I] = bestScore;
        curCol->PreviousVertex[I] = prevVertex;
        curCol->ReachingMove[I] = EndMove;
        return curCol;
    }


    const AlignmentColumn*
    PoaGraph::Impl::makeAlignmentColumn(Vertex v,
                                        const AlignmentColumnMap& alignmentColumnForVertex,
                                        const std::string& sequence,
                                        const PoaConfig& config)
    {
        AlignmentColumn* curCol = new AlignmentColumn(v, sequence.length() + 1);
        const PoaNode* vertexInfo = vertexInfoMap_[v];
        vector<const AlignmentColumn*> predecessorColumns =
                getPredecessorColumns(g_, v, alignmentColumnForVertex);

        //
        // handle read pos 0 separately:
        //
        if (predecessorColumns.size() == 0)
        {
            // if this vertex doesn't have any in-edges (^), then it has
            // no reaching move
            assert(v == enterVertex_);
            curCol->Score[0] = 0;
            curCol->ReachingMove[0] = InvalidMove;
            curCol->PreviousVertex[0] = null_vertex;
        }
        else if (config.UseLocalAlignment)
        {
            // under local alignment, we use the Start move
            curCol->Score[0] = 0;
            curCol->ReachingMove[0] = StartMove;
            curCol->PreviousVertex[0] = enterVertex_;
        }
        else
        {
            // otherwise it's a deletion
            float candidateScore;
            float bestScore = -FLT_MAX;
            Vertex prevVertex = null_vertex;
            MoveType reachingMove = InvalidMove;

            foreach (const AlignmentColumn * prevCol, predecessorColumns)
            {
                candidateScore = prevCol->Score[0] + config.Params.Missing;
                if (candidateScore > bestScore)
                {
                    bestScore = candidateScore;
                    prevVertex = prevCol->CurrentVertex;
                    reachingMove = DeleteMove;
                }
            }
            assert(reachingMove != InvalidMove);
            curCol->Score[0] = bestScore;
            curCol->ReachingMove[0] = reachingMove;
            curCol->PreviousVertex[0] = prevVertex;
        }

        //
        // tackle remainder of read.
        //
        // i represents position in array
        // readPos=i-1 represents position in read
        for (unsigned int i = 1, readPos = 0;  i <= sequence.length(); i++, readPos++)
        {
            float candidateScore;
            float bestScore = -FLT_MAX;
            Vertex prevVertex = null_vertex;
            MoveType reachingMove = InvalidMove;

            foreach (const AlignmentColumn* prevCol, predecessorColumns)
            {
                // Incorporate (Match or Mismatch)
                bool isMatch = sequence[readPos] == vertexInfo->Base;
                candidateScore = prevCol->Score[i - 1] + (isMatch ?
                                                             config.Params.Match :
                                                             config.Params.Mismatch);
                if (candidateScore > bestScore)
                {
                    bestScore = candidateScore;
                    prevVertex = prevCol->CurrentVertex;
                    reachingMove = (isMatch ? MatchMove : MismatchMove);
                }
                // Delete
                candidateScore = prevCol->Score[i] + config.Params.Missing;
                if (candidateScore > bestScore)
                {
                    bestScore = candidateScore;
                    prevVertex = prevCol->CurrentVertex;
                    reachingMove = DeleteMove;
                }
            }
            // Extra
            candidateScore = curCol->Score[i - 1] + config.Params.Extra;
            if (candidateScore > bestScore)
            {
                bestScore = candidateScore;
                prevVertex = v;
                reachingMove = ExtraMove;
            }
            assert(reachingMove != InvalidMove);
            curCol->Score[i] = bestScore;
            curCol->ReachingMove[i] = reachingMove;
            curCol->PreviousVertex[i] = prevVertex;
        }

        return curCol;
    }

    static void
    tagSpan(BoostGraph& g, Vertex start, Vertex end, VertexInfoMap& vertexInfoMap)
    {
        // cout << "Tagging span " << start << " to " << end << endl;
        std::list<Vertex> sortedVertices(num_vertices(g));
        topological_sort(g, sortedVertices.rbegin());
        bool spanning = false;
        foreach (Vertex v, sortedVertices)
        {
            if (v == start)
            {
                spanning = true;
            }
            if (v == end)
            {
                break;
            }
            if (spanning)
            {
                vertexInfoMap[v]->SpanningReads++;
            }
        }
    }

    static std::vector<Vertex>
    maxPath(BoostGraph& g, int totalReads, bool isLocal)
    {
        std::list<Vertex> path;
        VertexInfoMap vertexInfoMap = get(vertex_info, g);
        std::list<Vertex> sortedVertices(num_vertices(g));
        topological_sort(g, sortedVertices.rbegin());
        unordered_map<Vertex, Vertex> bestPrevVertex;

        // ignore ^ and $
        // TODO(dalexander): find a cleaner way to do this
        vertexInfoMap[sortedVertices.front()]->ReachingScore = 0;
        sortedVertices.pop_back();
        sortedVertices.pop_front();

        Vertex bestVertex = null_vertex;
        float bestReachingScore = -FLT_MAX;
        foreach (Vertex v, sortedVertices)
        {
            const PoaNode* vertexInfo = vertexInfoMap[v];
            int containingReads = vertexInfo->Reads;
            int spanningReads = vertexInfo->SpanningReads;
            float score = isLocal ?
                          (2 * containingReads - 1 * spanningReads - 0.0001f) :
                          (2 * containingReads - 1 * totalReads - 0.0001f);
            vertexInfoMap[v]->Score = score;
            vertexInfoMap[v]->ReachingScore = score;
            bestPrevVertex[v] = null_vertex;
            foreach (Edge e, in_edges(v, g))
            {
                Vertex sourceVertex = source(e, g);
                float rsc = score + vertexInfoMap[sourceVertex]->ReachingScore;
                if (rsc > vertexInfoMap[v]->ReachingScore)
                {
                    vertexInfoMap[v]->ReachingScore = rsc;
                    bestPrevVertex[v] = sourceVertex;
                }
                if (rsc > bestReachingScore)
                {
                    bestVertex = v;
                    bestReachingScore = rsc;
                }
            }
        }
        assert(bestVertex != null_vertex);

        // trace back from best-scoring vertex
        Vertex v = bestVertex;
        while (v != null_vertex)
        {
            path.push_front(v);
            v = bestPrevVertex[v];
        }
        return std::vector<Vertex>(path.begin(), path.end());
    }

    void PoaGraph::Impl::AddSequence(const std::string& sequence, const PoaConfig& config)
    {
        DEBUG_ONLY(repCheck());
        assert(sequence.length() > 0);

        const int I = sequence.length();
        int seqNo = sequences_.size();

        // Yes, it's true!  We need to retain a COPY of the SequenceFeatures
        // so that the shared_ptr's within will always have positive usage count
        // so long as this object exists.
        sequences_.push_back(sequence);

        if (seqNo == 0)
        {
            // first sequence in the alignment
            Vertex u = null_vertex, v;
            Vertex startSpanVertex = null_vertex, endSpanVertex;
            int readPos = 0;

            foreach (char base, sequence)
            {
                v = add_vertex(g_);
                vertexInfoMap_[v] = new PoaNode(base);
                if (readPos == 0)
                {
                    add_edge(enterVertex_, v, g_);
                    startSpanVertex = v;
                }
                else
                {
                    add_edge(u, v, g_);
                }
                u = v;
                readPos++;
            }
            assert(startSpanVertex != null_vertex);
            assert(u != null_vertex);
            endSpanVertex = u;
            add_edge(u, exitVertex_, g_);  // terminus -> $
            tagSpan(g_, startSpanVertex, endSpanVertex, vertexInfoMap_);
        }
        else
        {
            // calculate alignment column of sequence vs. graph
            AlignmentColumnMap alignmentColumnForVertex;
            vector<Vertex> sortedVertices(num_vertices(g_));
            topological_sort(g_, sortedVertices.rbegin());
            const AlignmentColumn* curCol;
            foreach (Vertex v, sortedVertices)
            {
                if (v != exitVertex_)
                {
                    curCol = makeAlignmentColumn(v, alignmentColumnForVertex,
                                                 sequence, config);
                }
                else
                {
                    curCol = makeAlignmentColumnForExit(v, alignmentColumnForVertex,
                                                        sequence, config);
                }
                alignmentColumnForVertex[v] = curCol;
            }

            // perform traceback from (I,$), threading the new sequence into the graph as
            // we go.
            int i = I;
            Vertex v = exitVertex_, forkVertex = exitVertex_;
            Vertex u = alignmentColumnForVertex[exitVertex_]->PreviousVertex[I];
            Vertex startSpanVertex, endSpanVertex = u;
            while ( !(u == enterVertex_ && i == 0) )
            {
                // u: current vertex
                // v: vertex last visited in traceback (could be == u)
                // forkVertex: the vertex that will be the target of a new edge
                int readPos = i - 1;
                curCol = alignmentColumnForVertex.at(u);
                assert(curCol != NULL);
                PoaNode* curNodeInfo = vertexInfoMap_[u];
                Vertex prevVertex = curCol->PreviousVertex[i];

                if (curCol->ReachingMove[i] == MatchMove)
                {
                    // if there is an extant forkVertex, join it
                    if (forkVertex != null_vertex)
                    {
                        add_edge(u, forkVertex, g_);
                        forkVertex = null_vertex;
                    }
                    // add to existing node
                    curNodeInfo->Reads++;
                    i--;
                }
                else if (curCol->ReachingMove[i] == DeleteMove ||
                         curCol->ReachingMove[i] == StartMove)
                {
                    if (forkVertex == null_vertex)
                    {
                        forkVertex = v;
                    }
                }
                else if (curCol->ReachingMove[i] == ExtraMove ||
                         curCol->ReachingMove[i] == MismatchMove)
                {
                    // begin a new arc with this read base
                    Vertex newForkVertex = add_vertex(g_);
                    vertexInfoMap_[newForkVertex] = new PoaNode(sequence[readPos]);
                    if (forkVertex != null_vertex)
                    {
                        add_edge(newForkVertex, forkVertex, g_);
                    }
                    else
                    {
                        add_edge(newForkVertex, v, g_);
                    }
                    forkVertex = newForkVertex;
                    i--;
                }
                else
                {
                    ShouldNotReachHere();
                }

                v = u;
                u = prevVertex;
            }
            startSpanVertex = v;
            if (startSpanVertex != exitVertex_)
            {
                tagSpan(g_, startSpanVertex, endSpanVertex, vertexInfoMap_);
            }

            // if there is an extant forkVertex, join it to enterVertex
            if (forkVertex != null_vertex)
            {
                add_edge(enterVertex_, forkVertex, g_);
                forkVertex = null_vertex;
            }

            // Clean up the mess we created.  Might be nicer to use scoped ptrs.
            foreach (AlignmentColumnMap::value_type& kv, alignmentColumnForVertex)
            {
                delete kv.second;
            }
        }

        DEBUG_ONLY(repCheck());
    }


    static boost::unordered_set<Vertex>
    childVertices(Vertex v,
                  BoostGraph& g)
    {
        boost::unordered_set<Vertex> result;
        foreach (Edge e, out_edges(v, g))
        {
            result.insert(target(e, g));
        }
        return result;
    }

    static boost::unordered_set<Vertex>
    parentVertices(Vertex v,
                   BoostGraph& g)
    {
        boost::unordered_set<Vertex> result;
        foreach (Edge e, in_edges(v, g))
        {
            result.insert(source(e, g));
        }
        return result;
    }


    tuple<string, float, vector< pair<Mutation*, float> >* >
    PoaGraph::Impl::FindConsensus(const PoaConfig& config)
    {
        std::stringstream ss;
        std::vector<Vertex> bestPath = maxPath(g_, NumSequences(), config.UseLocalAlignment);
        foreach (Vertex v, bestPath)
        {
            PoaNode* consensusNode = vertexInfoMap_[v];
            consensusNode->IsInConsensus = true;
            ss << consensusNode->Base;
        }

        // if requested, identify likely sequence variants

        // will be deallocated by PoaConsensus destructor.
        vector< pair<Mutation*, float> >* variants = new vector< pair<Mutation*, float> >();

        if (true)  // TODO(dalexander): Add a flag to PoaConfig
        {
            for (int i = 2; i < (int)bestPath.size() - 2; i++) // NOLINT
            {
                Vertex v = bestPath[i];
                boost::unordered_set<Vertex> children = childVertices(v, g_);

                // Look for a direct edge from the current node to the node
                // two spaces down---suggesting a deletion with respect to
                // the consensus sequence.
                if (children.find(bestPath[i + 2]) != children.end())
                {
                    float score = -vertexInfoMap_[bestPath[i + 1]]->Score;
                    variants->push_back(make_pair(new Mutation(DELETION, i + 1, '-'), score));
                }

                // Look for a child node that connects immediately back to i + 1.
                // This indicates we should try inserting the base at i + 1.

                // Parents of (i + 1)
                boost::unordered_set<Vertex> lookBack = parentVertices(bestPath[i + 1], g_);

                // (We could do this in STL using std::set sorted on score, which would then
                // provide an intersection mechanism (in <algorithm>) but that actually ends
                // up being more code.  Sad.)
                float bestInsertScore = -FLT_MAX;
                Vertex bestInsertVertex = null_vertex;

                foreach (Vertex v, children)
                {
                    boost::unordered_set<Vertex>::iterator found = lookBack.find(v);
                    if (found != lookBack.end())
                    {
                        float score = vertexInfoMap_[*found]->Score;
                        if (score > bestInsertScore)
                        {
                            bestInsertScore = score;
                            bestInsertVertex = *found;
                        }
                    }
                }

                if (bestInsertVertex != null_vertex)
                {
                    char base = vertexInfoMap_[bestInsertVertex]->Base;
                    variants->push_back(
                            make_pair(new Mutation(INSERTION, i + 1, base), bestInsertScore));
                }

                // Look for a child node not in the consensus that connects immediately
                // to i + 2.  This indicates we should try mismatching the base i + 1.

                // Parents of (i + 2)
                lookBack = parentVertices(bestPath[i + 2], g_);

                float bestMismatchScore = -FLT_MAX;
                Vertex bestMismatchVertex = null_vertex;

                foreach (Vertex v, children)
                {
                    if (v == bestPath[i + 1]) continue;

                    boost::unordered_set<Vertex>::iterator found = lookBack.find(v);
                    if (found != lookBack.end())
                    {
                        float score = vertexInfoMap_[*found]->Score;
                        if (score > bestMismatchScore)
                        {
                            bestMismatchScore = score;
                            bestMismatchVertex = *found;
                        }
                    }
                }

                if (bestMismatchVertex != null_vertex)
                {
                    // TODO(dalexander): As implemented (compatibility), this returns
                    // the score of the mismatch node. I think it should return the score
                    // difference, no?
                    char base = vertexInfoMap_[bestMismatchVertex]->Base;
                    variants->push_back(
                            make_pair(new Mutation(SUBSTITUTION, i + 1, base), bestMismatchScore));
                }
            }
        }

        return boost::make_tuple(ss.str(), 0.0f, variants);  // TODO(dalexander): where do we get scores?
    }

    inline int
    PoaGraph::Impl::NumSequences() const
    {
        return sequences_.size();
    }

    string PoaGraph::Impl::ToGraphViz(int flags) const
    {
        std::stringstream ss;
        write_graphviz(ss, g_, my_label_writer(vertexInfoMap_,
                                               flags & COLOR_NODES,
                                               flags & VERBOSE_NODES));
        return ss.str();
    }


    void
    PoaGraph::Impl::WriteGraphVizFile(string filename, int flags) const
    {
        std::ofstream outfile(filename.c_str());
        outfile << ToGraphViz(flags);
        outfile.close();
    }

    // PIMPL idiom delegation

    void
    PoaGraph::AddSequence(const std::string& sequence, const PoaConfig& config)
    {
        impl->AddSequence(sequence, config);
    }

    int
    PoaGraph::NumSequences() const
    {
        return impl->NumSequences();
    }

    tuple<string, float, std::vector< std::pair<Mutation*, float> >* >
    PoaGraph::FindConsensus(const PoaConfig& config) const
    {
        return impl->FindConsensus(config);
    }

    string
    PoaGraph::ToGraphViz(int flags) const
    {
        return impl->ToGraphViz(flags);
    }

    void
    PoaGraph::WriteGraphVizFile(string filename, int flags) const
    {
        impl->WriteGraphVizFile(filename, flags);
    }

    PoaGraph::PoaGraph()
    {
        impl = new Impl();
    }

    PoaGraph::~PoaGraph()
    {
        delete impl;
    }
}


