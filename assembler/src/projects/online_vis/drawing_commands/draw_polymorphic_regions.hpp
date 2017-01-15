//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../environment.hpp"
#include "../command.hpp"
#include "../errors.hpp"
#include "io/reads/wrapper_collection.hpp"

namespace online_visualization {

class DrawPolymorphicRegions : public DrawingCommand {

    GraphComponent<Graph> ConstructComponent(DebruijnEnvironment& curr_env, EdgeId firstEdge, EdgeId secondEdge, size_t windowSize) const
    {
        PathStorageCallback<Graph> callback(curr_env.graph());
        ProcessPaths(curr_env.graph(), 0, windowSize*4,
                        curr_env.graph().EdgeEnd(firstEdge), curr_env.graph().EdgeStart(secondEdge),
                        callback);
        vector<vector<EdgeId>> paths = callback.paths();
        vector<VertexId> verticesToAdd;
        verticesToAdd.push_back(curr_env.graph().EdgeEnd(firstEdge));
        for(auto edges : paths)
        {
            for(auto edge : edges)
            {
                verticesToAdd.push_back(curr_env.graph().EdgeEnd(edge));
            }
        }
        return GraphComponent<Graph>::FromVertices(curr_env.graph(), verticesToAdd);
    }

    void DrawPicture(DebruijnEnvironment& curr_env, Sequence& genome) const {
        size_t windowSize = 500;
        for(size_t i = 0; i < genome.size() - windowSize - 1 - curr_env.k_value(); ++i)
        {
            RtSeq firstKmer = genome.Subseq(i).start<RtSeq>(curr_env.k_value() + 1);
            RtSeq secondKmer = genome.Subseq(i + windowSize).start<RtSeq>(curr_env.k_value() + 1);
            firstKmer = curr_env.kmer_mapper().Substitute(firstKmer);
            secondKmer = curr_env.kmer_mapper().Substitute(secondKmer);
            pair<EdgeId, size_t> positionFirst = curr_env.index().get(firstKmer);
            if(positionFirst.first == EdgeId(0))
            {
                continue;
            }

            if(curr_env.graph().length(positionFirst.first) < 300)
            {
                i += curr_env.graph().length(positionFirst.first) - positionFirst.second;
                continue;
            }
            else
            {
                pair<EdgeId, size_t> positionSecond = curr_env.index().get(secondKmer);
                if(positionSecond.first == EdgeId(0))
                {
                    continue;
                }

                if(curr_env.graph().length(positionSecond.first) < 300)
                {
                    i += curr_env.graph().length(positionSecond.first) - positionSecond.second;
                    continue;
                }
                else
                {
                    if(positionFirst.first == positionSecond.first)
                    {
                        i += curr_env.graph().length(positionSecond.first) - positionSecond.second;
                        continue;
                    }
                    INFO("Constructing component around " << i << "-th position in the genome");
                    GraphComponent<Graph> polymorphicRegion = ConstructComponent(curr_env, positionFirst.first, positionSecond.first, windowSize);

                    if(polymorphicRegion.e_size() > 5)
                    {
                        using namespace visualization::visualization_utils;
                        WriteComponentSinksSources(polymorphicRegion,
                                                   curr_env.folder() + "/" +
                                                           ToString(curr_env.graph().int_id(*polymorphicRegion.vertices().begin())) + ".dot",
                                                   visualization::graph_colorer::DefaultColorer(curr_env.graph()),
                                                   curr_env.labeler());

                        INFO("Component is written to " + curr_env.folder() + ToString(curr_env.graph().int_id(*polymorphicRegion.vertices().begin())) + ".dot");
                    }

                    i += curr_env.graph().length(positionSecond.first) - positionSecond.second;
                    continue;
                }
            }

        }
    }


protected:
    size_t MinArgNumber() const {
        return 0;
    }

    bool CheckCorrectness(const vector<string>&) const {
        return true;
    }

public:
    string Usage() const {
        string answer;
        answer = answer + "Command `draw_polymorphic` \n" +
                        "Usage:\n" +
                        "> draw_polymorphic\n" +
                        " You should run load_genome command before it, to proceed. \n" +
                        "This command draws polymorphic regions between two conserved edges.";
        return answer;
    }

    DrawPolymorphicRegions() : DrawingCommand("draw_polymorphic")
    {
    }

    void Execute(DebruijnEnvironment& curr_env, const ArgumentList& arg_list) const {
        const vector<string>& args = arg_list.GetAllArguments();
        if (!CheckCorrectness(args))
            return;
        make_dir(curr_env.folder());
        Sequence genome = curr_env.genome();
        if(genome.size() == 0)
        {
            cout << "Reference genome is not uploaded\n";
            return;
        }
        DrawPicture(curr_env, genome);
        INFO("End");

    }
};
}

