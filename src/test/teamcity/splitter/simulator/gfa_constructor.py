#!/usr/bin/env python3
import argparse
import itertools
import shutil
import random
import os

from collections import namedtuple

Segment = namedtuple('Segment', ['id', 'length', 'seq'])
Vertex = namedtuple('Vertex', ['id', 'overlap'])
Link = namedtuple('Link', ['start', 'end', 'overlap'])
Path = namedtuple('Path', ['id', 'edges'])


NUCL_LIST = ['A', 'C', 'G', 'T']
RC = {'A': 'T',
      'C': 'G',
      'G': 'C',
      'T': 'A'}

def generate_random_sequence(length):
    return "".join([random.choice(NUCL_LIST) for _ in range(length)])


class OverlapGraph(object):
	def __init__(self):
		self.edges = {}
		self.vertices = {}
		#edge to end-overlap
		self.ends = {}
		#edge to start-overlap
		self.starts = {}
		self.vertex_to_out = {}
		self.vertex_to_in = {}
		self.links = []
		self.paths = []
		# self.edge_id = 0
		# self.vertex_id = 0


	def add_segment(self, edge_id, length):
		segment = Segment(id=edge_id, length=length, seq="")
		self.edges[edge_id] = segment
		# edge_id += 1


	def add_vertex(self, vertex_id, overlap):
		vertex = Vertex(id=vertex_id, overlap=overlap)
		self.vertices[vertex_id] = vertex


	def add_link(self, first_id, second_id):
		if first_id in self.vertices and second_id in self.edges:
			self._add_outlink(vertex_id=first_id, edge_id=second_id)
		elif first_id in self.edges and second_id in self.vertices:
			self._add_inlink(vertex_id=second_id, edge_id=first_id)
		else:
			raise ValueError("Invalid id pair {} {}".format(first_id, second_id))


	def add_gfa_link(self, first_id, second_id, overlap):
		link = Link(start=first_id, end=second_id, overlap=overlap)
		self.links.append(link)


	#simple graph: every vertex has the same in- and out- degree
	def check_simplicity(self):
		for vertex in self.vertices.keys():
			assert(len(self.vertex_to_in[vertex]) == len(self.vertex_to_out[vertex]))
		for edge, end in self.ends.items():
			if edge in self.starts:
				start = self.starts[edge]
				in_overlap = self.vertices[end].overlap
				out_overlap = self.vertices[start].overlap
				edge_len = self.edges[edge].length
				assert(edge_len >= in_overlap + out_overlap)


	#check if link set corresponds to the graph structure (end of first edge == start of second edge)
	def check_link_vertex_consistency(self):
		for link in self.links:
			print(self.starts)
			end_first = self.ends[link.start]
			start_second = self.starts[link.end]
			assert(start_second == end_first)
			assert(link.overlap == self.vertices[start_second].overlap)


	def construct_sequences(self):
		vertex_to_seq = {}
		for edge in self.edges.keys():
			if edge in self.ends:
				end_vertex = self.ends[edge]
				if end_vertex not in vertex_to_seq:
					vertex_len = self.vertices[end_vertex].overlap
					vertex_seq = generate_random_sequence(vertex_len)
					vertex_to_seq[end_vertex] = vertex_seq
			if edge in self.starts:
				start_vertex = self.starts[edge]
				if start_vertex not in vertex_to_seq:
					vertex_len = self.vertices[start_vertex].overlap
					vertex_seq = generate_random_sequence(vertex_len)
					vertex_to_seq[start_vertex] = vertex_seq

		new_edges = {}
		for edge in self.edges:
			start_seq = ""
			end_seq = ""
			if edge in self.starts:
				start_seq = vertex_to_seq[self.starts[edge]]
			if edge in self.ends:
				end_seq = vertex_to_seq[self.ends[edge]]
			mid_seq_len = self.edges[edge].length - len(start_seq) - len(end_seq)
			assert(mid_seq_len >= 0)
			mid_seq = generate_random_sequence(mid_seq_len)
			edge_seq = start_seq + mid_seq + end_seq
			new_edges[edge] = Segment(id=edge, length=self.edges[edge].length, seq=edge_seq)
			# print(edge.seq)
		self.edges = new_edges


	def add_path(self, path):
		self.paths.append(path)


	def check_paths(self):
		print(self.starts)
		print(self.ends)
		covered_edges = set()
		for path in self.paths:
			for edge in path.edges:
				print(edge)
				if edge in covered_edges:
					raise ValueError("Edge {} covered by multiple paths".format(edge))
				covered_edges.add(edge)
		for path in self.paths:
			for first_edge, second_edge in itertools.pairwise(path.edges):
				if self.ends[first_edge] != self.starts[second_edge]:
					raise ValueError("Edges {} and {} are not incident".format(first_edge, second_edge))
		return True


	#link from end of the segment to the vertex
	def _add_inlink(self, vertex_id, edge_id):
		assert(edge_id not in self.ends)
		self.ends[edge_id] = vertex_id
		if vertex_id not in self.vertex_to_in:
			self.vertex_to_in[vertex_id] = []
		self.vertex_to_in[vertex_id].append(edge_id)


	#Link from the vertex to the start of the segment
	def _add_outlink(self, vertex_id, edge_id):
		assert(edge_id not in self.starts)
		self.starts[edge_id] = vertex_id
		if vertex_id not in self.vertex_to_out:
			self.vertex_to_out[vertex_id] = []
		self.vertex_to_out[vertex_id].append(edge_id)


#Graph format descripition
# Graph name
# Number of segments/edges
# List of segments <id, length>
# Number of overlaps/vertices
# List of overlaps <id, length>
# Number of edge-vertex links
# List of edge-vertex links <vertexid, segmentid> or <segmentid, linkid>
# Number of edge-edge links
# List of edge-edge links <edgeid, edgeid>
# Number of paths
# List of paths <sequence of edges>


def parse_graph(input_path):
	graph = OverlapGraph()
	with open(input_path, "r") as in_handle:
		graph_name = in_handle.readline().strip()
		num_segments = int(in_handle.readline().strip())
		for _ in range(num_segments):
			line = in_handle.readline().strip().split()
			edge_id, length = line[0], int(line[1])
			graph.add_segment(edge_id, length)
		num_overlaps = int(in_handle.readline())
		for _ in range(num_overlaps):
			line = in_handle.readline().strip().split()
			vertex_id, length = line[0], int(line[1])
			graph.add_vertex(vertex_id, length)
		num_links = int(in_handle.readline())
		for _ in range(num_links):
			line = in_handle.readline().strip().split()
			first_id, second_id = line[0], line[1]
			graph.add_link(first_id, second_id)
		num_gfa_links = int(in_handle.readline())
		for _ in range(num_gfa_links):
			line = in_handle.readline().strip().split()
			first_id, second_id, overlap = line[0], line[1], int(line[2])
			graph.add_gfa_link(first_id, second_id, overlap)
		num_paths = int(in_handle.readline())
		for _ in range(num_paths):
			line = in_handle.readline().strip().split()
			print(line)
			name = line[0]
			edges = line[1:]
			print(edges)
			path = Path(id=name, edges=edges)
			graph.add_path(path)
	return graph_name, graph


def output_graph(graph, name, output_path):
	gfa_path = os.path.join(output_path, name + ".gfa")
	ctg_path = os.path.join(output_path, name + ".fasta")
	default_cov = 1.0
	with open(gfa_path, "w") as gfa_handle, open(ctg_path, "w") as ctg_handle:
		gfa_handle.write("H\tVN:Z:1.0\n")
		for edge in graph.edges.values():
			print(edge)
			gfa_handle.write("S\t{}\t{}\n".format(edge.id, edge.seq))
		for link in graph.links:
			gfa_handle.write("L\t{}\t+\t{}\t+\t{}M\n".format(link.start, link.end, link.overlap))
		# for vertex in graph.vertices.keys():
		# 	if vertex in graph.vertex_to_in and vertex in graph.vertex_to_out:
		# 		vertex_len = graph.vertices[vertex].overlap
		# 		for in_edge in graph.vertex_to_in[vertex]:
		# 			for out_edge in graph.vertex_to_out[vertex]:
		# 				gfa_handle.write("L\t{}\t+\t{}\t+\t{}M\n".format(in_edge, out_edge, vertex_len))
		for path in graph.paths:
			edge_string = "+,".join(path.edges) + "+"
			ov_lengths = []
			for edge in path.edges[:-1]:
				end_vertex = graph.ends[edge]
				vertex_len = graph.vertices[end_vertex].overlap
				ov_lengths.append(vertex_len)
			vertex_string = ",".join(["{}M".format(length) for length in ov_lengths])
			gfa_handle.write("P\t{}\t{}\t{}\n".format(path.id, edge_string, vertex_string))

		for path in graph.paths:
			ctg_seq = graph.edges[path.edges[0]].seq
			for edge in path.edges[1:]:
				start_vertex = graph.starts[edge]
				vertex_len = graph.vertices[start_vertex].overlap
				edge_seq = graph.edges[edge].seq
				ctg_seq += edge_seq[vertex_len:]
			ctg_handle.write(">{}\n{}\n".format(path.id, ctg_seq))


def createparser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help="Graph description in text form (num segments, list of segments, num overlaps, list of overlaps, list of links)")
	parser.add_argument('--output', help="Output directory")
	return parser


if __name__ == "__main__":
	parser = createparser()
	args = parser.parse_args()

	output_path = args.output

	if os.path.exists(output_path):
		shutil.rmtree(output_path)
	os.mkdir(output_path)

	name, graph = parse_graph(args.input)
	print("Checking simplicity")
	graph.check_simplicity()
	graph.check_link_vertex_consistency()
	graph.construct_sequences()
	print("Checking paths")
	graph.check_paths()
	print("Outputting graph")
	output_graph(graph, name, args.output)