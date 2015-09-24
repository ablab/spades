//
// Created by andrey on 21.09.15.
//

#include "scaffold_graph_visualizer.hpp"

namespace scaffold_graph {

const map<size_t, string> ScaffoldEdgeColorer::color_map =
        {{(size_t) -1, "black"},
         {0, "red"},
         {1, "blue"},
         {2, "green"},
         {3, "magenta"},
         {4, "orange"},
         {5, "cyan"}};

const string ScaffoldEdgeColorer::default_color = "black";

}