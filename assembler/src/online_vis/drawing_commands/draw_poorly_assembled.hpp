//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "sequence_mapper.hpp"

namespace online_visualization {

class PathNeighbourhoodFinder : public AbstractNeighbourhoodFinder<Graph> {
private:

    VertexId OtherEnd(EdgeId e, VertexId v) const {
        if (this->graph().EdgeStart(e) == v)
            return this->graph().EdgeEnd(e);
        else
            return this->graph().EdgeStart(e);
    }

    bool Go(VertexId v, size_t curr_depth, set<VertexId>& grey, set<VertexId>& black) const {
        //allows single vertex be visited many times with different depth values
        if (curr_depth >= max_depth_)
            return true;
        if (grey.size() >= max_size_)
            return false;

        grey.insert(v);
        vector<EdgeId> incident_path;
        vector<EdgeId> incident_non_path;
        for (EdgeId e : graph().IncidentEdges(v)) {
            if (path_edges_.count(e) && /*condition to stretch forward*/ graph().EdgeStart(e) == v)
                incident_path.push_back(e);
            else
                incident_non_path.push_back(e);
        }

        for (EdgeId e : incident_non_path) {
            if (graph().length(e) > edge_length_bound_)
                continue;
            if (!Go(OtherEnd(e, v), curr_depth + 1, grey, black))
                return false;
        }

        for (EdgeId e : incident_path) {
            if (!Go(OtherEnd(e, v), curr_depth, grey, black))
                return false;
        }

        black.insert(v);
        return true;
    }

public:
    static const size_t DEFAULT_EDGE_LENGTH_BOUND = 500;
    static const size_t DEFAULT_MAX_DEPTH = 3;
    static const size_t DEFAULT_MAX_SIZE = 30;

    set<EdgeId> path_edges_;
    const size_t edge_length_bound_;
    const size_t max_depth_;
    const size_t max_size_;

    set<VertexId> last_inner_;

    PathNeighbourhoodFinder(const Graph &graph, const vector<EdgeId>& path, size_t edge_length_bound = DEFAULT_EDGE_LENGTH_BOUND,
                            size_t max_depth = DEFAULT_MAX_DEPTH, size_t max_size = DEFAULT_MAX_SIZE)
            : AbstractNeighbourhoodFinder<Graph>(graph),
              path_edges_(path.begin(), path.end()),
              edge_length_bound_(edge_length_bound),
              max_depth_(max_depth),
              max_size_(max_size) {
    }


    GraphComponent<Graph> Find(VertexId v) {
        last_inner_.clear();
        set<VertexId> grey;
        set<VertexId> black;
        Go(v, 0, grey, black);
        last_inner_ = black;
        last_inner_.insert(v);
        return GraphComponent<Graph>(graph(), grey.begin(), grey.end());
    }

    vector<VertexId> InnerVertices(const GraphComponent<Graph> &/*component*/) {
        return vector<VertexId>(last_inner_.begin(), last_inner_.end());
    }
};

class DrawPoorlyAssembledCommand : public DrawingCommand {
    const double WELL_ASSEMBLED_CONSTANT = 0.7;
private:

    shared_ptr<GraphSplitter<Graph>> SplitterAlongPath(
            const Graph &graph, const Path<EdgeId>& path) const {
        typedef typename Graph::VertexId VertexId;
        shared_ptr<RelaxingIterator<VertexId>> inner_iterator = make_shared<
                PathIterator<Graph>>(graph, path);
        shared_ptr<AbstractNeighbourhoodFinder<Graph>> nf = make_shared<PathNeighbourhoodFinder>(graph, path.sequence());
        return make_shared<NeighbourhoodFindingSplitter<Graph>>(graph,
                                                                inner_iterator, nf);
    }

    void WriteAlongPath(const Graph& g, Path<EdgeId> path, const string& folder_name,
                                  shared_ptr<omnigraph::visualization::GraphColorer<Graph>> colorer,
                                  const GraphLabeler<Graph> &labeler, bool color_path = true) const {
        using namespace omnigraph::visualization;
        using std::make_shared;
        auto edge_colorer = make_shared<CompositeEdgeColorer<Graph>>("black");
        edge_colorer->AddColorer(colorer);
        if (color_path) {
            edge_colorer->AddColorer(make_shared<SetColorer<Graph>>(g, path.sequence(), "green"));
        }
        shared_ptr<GraphColorer<Graph>> resulting_colorer = make_shared<CompositeGraphColorer<Graph>>(colorer, edge_colorer);
        shared_ptr<GraphSplitter<Graph>> rs = SplitterAlongPath(g, path);
        auto filter = make_shared<omnigraph::SmallComponentFilter<Graph>>(g, 3);
        shared_ptr<GraphSplitter<Graph>> splitter = make_shared<omnigraph::CondensingSplitterWrapper<Graph>>(rs, filter);
        WriteComponents<Graph>(g, folder_name, splitter, resulting_colorer, labeler);
    }

    void DrawPicturesAlongGenomePart(DebruijnEnvironment& curr_env, const Sequence& piece_of_genome, string label = "") const {
        const MappingPath<EdgeId>& mapping_path = curr_env.mapper().MapSequence(piece_of_genome);
        make_dir(curr_env.folder());
        stringstream namestream;
        namestream << curr_env.folder() << "/" << curr_env.GetFormattedPictureCounter()
                << "_" << curr_env.file_name() << "/";
        make_dir(namestream.str());
        namestream << label;
        make_dir(namestream.str());
        WriteAlongPath(curr_env.graph(), mapping_path.path(), namestream.str(),
                                  curr_env.coloring(), curr_env.labeler());

        cout << "The pictures is written to " << namestream.str() << endl;

        curr_env.inc_pic_counter();
    }

    void DrawContig(DebruijnEnvironment& curr_env, io::SingleRead contig) const {
        Sequence seq = contig.sequence();
        string label = contig.name();
        DrawPicturesAlongGenomePart(curr_env, seq, label);
        LOG("Contig " << contig.name() << " has been drawn");
    }

    io::SingleRead MakeValid(const io::SingleRead& contig) const {
        std::string str = contig.GetSequenceString();
        for (size_t i = 0; i < str.length(); ++i) {
            if (str[i] == 'N')
                str[i] = nucl(char(i % 4));
        }
        return io::SingleRead(contig.name(), str);
    }

    bool IsPoorlyAssembled(const GraphPack& gp, io::SingleRead contig, string base_assembly_prefix) const {
        MappingPath<EdgeId> mapping_path = debruijn_graph::MapperInstance(gp)->MapRead(contig);
        auto pos_handler = gp.edge_pos;
        map<string, size_t> base_ctg_2_len;
        for (EdgeId e : mapping_path.simple_path()) {
            auto positions = pos_handler.GetEdgePositions(e);
            for (EdgePosition pos : positions) {
                if (boost::starts_with(pos.contigId, base_assembly_prefix)) {
                    base_ctg_2_len[pos.contigId] += pos.mr.mapped_range.size();
                }
            }
        }
        for (pair<string, size_t> entry : base_ctg_2_len)
            if (double(entry.second) > double(contig.size()) * WELL_ASSEMBLED_CONSTANT)
                return false;

        return true;
    }

protected:
    size_t MinArgNumber() const {
        return 2;
    }

    bool CheckCorrectness(const vector<string>& args) const {
        if (!CheckEnoughArguments(args))
            return false;

        return true;
    }

public:
    string Usage() const {
        string answer;
        answer = answer + "Command `draw_poorly_assembled` \n" + "Usage:\n"
                + "> draw_poorly_assembled <contigs_file> <prefix_of_base_assembly> [first N contigs to analyze]\n"
                + " Draws pictures of contigs that are not well covered with any contig in base assembly.";
        return answer;
    }

    DrawPoorlyAssembledCommand()
            : DrawingCommand("draw_poorly_assembled") {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;

        std::string contigs_file = args[1];
        string base_assembly_prefix = args[2];

        if (!CheckFileExists(contigs_file)) {
            LOG("File with contigs " << contigs_file << " not found");
        }

        size_t contig_cnt = -1u;
        if (args.size() > 3) {
            LOG("Will analyze first " << args[3] << " contigs");
            contig_cnt = lexical_cast<size_t>(args[3]);
        }

        io::FileReadStream irs(contigs_file);

        size_t i = 0;
        while (!irs.eof() && i < contig_cnt) {
            io::SingleRead contig;
            irs >> contig;
            contig = MakeValid(contig);
            LOG("Considering contig " << contig.name());

            // if read is valid and also the name contains a given string <contig_name> as a substring.
            VERIFY(contig.IsValid());

            if (IsPoorlyAssembled(curr_env.graph_pack(), contig, base_assembly_prefix)) {
                LOG("Was poorly assembled, drawing");
                DrawContig(curr_env, contig);
            } else {
                LOG("Was well assembled");
            }

            ++i;
        }

    }

};
}
