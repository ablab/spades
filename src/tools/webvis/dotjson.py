
############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import dot_parser
import jsonpickle
from random import randint

class JNode(object):
    def __init__(self, id, label=None):
        if id is None:
            print "Unnamed node detected!"
            self.id = "[unknown]"
        self.id = id
        if label is None:
            self.label = id
        else:
            self.label = label
        self.x = randint(0, 10)
        self.y = randint(0, 10)
        self.size = 3

    @staticmethod
    def from_node(node):
        return JNode(
            node.get_name(),
            #node.get_label()
        )

    def to_json(self):
        return self.__dict__

class JEdge(object):
    def __init__(self, source, target, id=None, label=None):
        if id is None:
            self.id = source + "->" + target
        else:
            self.id = id
        self.source = source
        self.target = target
        if label is not None:
            self.label = label

    @staticmethod
    def from_edge(edge):
        return JEdge(
            edge.get_source().split(":")[0],
            edge.get_destination().split(":")[0],
            edge.get_id(),
            edge.get_label()
        )

class JGraph(object):

    def __init__(self, nodes, edges):
        self.nodes = nodes
        unique_edges = dict()
        for edge in edges:
            unique_edges[edge.id] = edge
        self.edges = unique_edges.values()
    
    @staticmethod
    def from_graph(graph):
        return JGraph(
            map(JNode.from_node, graph.get_nodes()[1:]),
            map(JEdge.from_edge, graph.get_edges())
        )
    
def dot_to_json(graph):
    graph = JGraph.from_graph(dot_parser.parse_dot_data(graph))
    return jsonpickle.encode(graph)
