############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from abstract_graph import Abstract_Graph, Abstract_Edge, Abstract_Vertex
import abstract_graph
import os
import saveparser
import experimental
import graph
import logging
import cStringIO
import itertools
from test_util import TestUtils
from rectangle import Rectangle
from utils import conjugate
from collections import defaultdict
import utils

DEEP_THRESHOLD = 10

def avoid_N(x, y):
    if x != 'N':
        return x
    return y


class BVertex(Abstract_Vertex):
    vid = 0

    def __init__(self, key):
        Abstract_Vertex.__init__(self, BVertex.vid)
        self.key = key
        BVertex.vid += 1


class BEdge(Abstract_Edge):
    eid = 0

    def __init__(self, v1, v2, diag):
        Abstract_Edge.__init__(self, v1, v2, BEdge.eid)
        self.diagonals = [diag]
        BEdge.eid += 1

    def length(self): 
        return reduce(lambda l, d: l + d.offsetc - d.offseta, self.diagonals, 0)
        
    
    def _get_splitted_seq(self, K, d, is_sc): 
        (seq1, seq2) = self.get_paired_seq(K, d)
        seq = ''.join(map(avoid_N, seq1, seq2)).strip('N')
        return  seq.split(self.get_midle_seq())
    
    def get_begin_seq(self, K, d, is_sc):
        return self._get_splitted_seq(K, d, is_sc)[0]

    def get_end_seq(self, K, d, is_sc):
        return self._get_splitted_seq(K, d, is_sc)[1] 

    def get_midle_seq(self):
        return reduce(lambda seq, d: seq + d.rectangle.e1.seq[d.offseta:d.offsetc], self.diagonals, "")
    
    def __can_add(self, d, step, init_index, CUT_LENGTH_THRESHOLD, CUT_THRESHOLD):
        cur_len = 0
        diag_index = init_index
        diag = self.diagonals[diag_index]
        while cur_len < d:
            if diag.offsetc - diag.offseta < CUT_LENGTH_THRESHOLD or (diag.support() < CUT_THRESHOLD):
                return False   
            diag_index += step
            if diag_index == -1 or diag_index == len(self.diagonals):
                return True
            diag = self.diagonals[diag_index]
            cur_len += diag.offsetc - diag.offseta
        return True

    def get_seq_for_contig(self, K, d, is_sc, is_careful):
        if is_sc:
            CUT_THRESHOLD = 2.0 #TODO: should take from histogram
            CUT_LENGTH_THRESHOLD = 5.0
            MIN_LENGTH = 4 * d
        else:
            CUT_LENGTH_THRESHOLD = 5
            CUT_THRESHOLD = 0.0
            MIN_LENGTH = 0
        (seq1, seq2) = self.get_paired_seq(K, d)
        seq = ''.join(map(avoid_N, seq1, seq2)).strip('N')
        first = self.diagonals[0]
        last = self.diagonals[-1]
        if len(seq1) > MIN_LENGTH:
            can_add_begin = False;
            can_add_end = False;
            if (not is_careful):
                can_add_begin = self.__can_add(d, 1, 0, CUT_LENGTH_THRESHOLD, CUT_THRESHOLD)
                can_add_end = self.__can_add(d, -1, len(self.diagonals) - 1, CUT_LENGTH_THRESHOLD, CUT_THRESHOLD)
            begin = first.rectangle.e1.seq[:first.offseta] if can_add_begin else ''
            end = last.rectangle.e2.seq[last.offsetd + K:] if can_add_end else ''
            if can_add_end or can_add_begin:
                return begin + seq + end

        seq1 = cStringIO.StringIO()
        for this in self.diagonals:
            seq1.write(this.rectangle.e1.seq[this.offseta: this.offsetc])
        last = self.diagonals[-1]
        seq1.write(last.rectangle.e1.seq[last.offsetc:])#last.offsetc + d])
        first = self.diagonals[0]
        seq1 = first.rectangle.e2.seq[:first.offsetb] + seq1.getvalue()[
                                                        d:]#[first.offsetb - d:first.offsetb] + seq1.getvalue()[d:]
        return seq1

    def get_seq(self, K, d):
        (seq1, seq2) = self.get_paired_seq(K, d)
        seq = utils.seq_join(seq1, seq2).strip('N')
        return seq

    def get_paired_seq(self, K, d):
        seq1 = cStringIO.StringIO()
        seq2 = cStringIO.StringIO()
        seq2.write('N' * d)
        for this in self.diagonals:
            seq1.write(this.rectangle.e1.seq[this.offseta: this.offsetc])
            seq2.write(this.rectangle.e2.seq[this.offsetb: this.offsetd])
        last = self.diagonals[-1]
        seq1.write(last.rectangle.e1.seq[last.offsetc: last.offsetc + K])
        seq2.write(last.rectangle.e2.seq[last.offsetd: last.offsetd + K])
        seq1.write('N' * (len(seq2.getvalue()) - len(seq1.getvalue())))
        return (seq1.getvalue(), seq2.getvalue())

    def get_cvr(self):
        sumlen = reduce(lambda s, d: s + d.offsetc - d.offseta, self.diagonals, 0)
        cvr = reduce(lambda cvr, d: (d.rectangle.e1.cvr + d.rectangle.e2.cvr) * 0.5 * (d.offsetc - d.offseta), self.diagonals, 0.0) 
        cvr /= sumlen
        return cvr

    def __repr__(self):
        return str((self.eid, self.diagonals))


class BGraph(Abstract_Graph):
    def __init__(self, graph, d, test_utils):
        Abstract_Graph.__init__(self)
        self.logger = logging.getLogger('rectangles')
        self.graph = graph
        self.d = d
        self.diagonals = set()
        self.test_utils = test_utils

    def __remove_bedge__(self, bedge):
        bv1 = bedge.v1
        bv2 = bedge.v2
        if bedge in bv1.out:
            bv1.out.remove(bedge)
        if (bedge in bv2.inn):
            bv2.inn.remove(bedge)
        self.__try_delete_bv(bv1)
        self.__try_delete_bv(bv2)
        if bedge.eid in self.es:
            del self.es[bedge.eid]
        for diag in bedge.diagonals:
            if diag in self.diagonals:
                self.diagonals.remove(diag)
    
    def __remove_bedge_and_conj(self, bedge):
        self.__remove_bedge__(bedge)
        self.__remove_bedge__(bedge.conj)

    def use_scaffold_paired_info(self, L, additional_prd):
        long_edges = set()
        used_paires = set()
        connect_edges = set()
        count_correct_scaffolds = 0
        count_incorrect_scaffolds = 0
        for edge_id, edge in self.es.items():
            if edge.length() > L:
                long_edges.add(edge)
        for e1 in long_edges:
            for e2 in long_edges:
                first_rectangle = e1.diagonals[-1].rectangle
                second_rectangle = e2.diagonals[0].rectangle
                e11 = first_rectangle.e1
                e12 = first_rectangle.e2
                e21 = second_rectangle.e1
                e22 = second_rectangle.e2
                if (e12.eid,
                    e21.eid) in additional_prd: #or (e11, e22) in additional_prd or (e11, e21) in additional_prd or (e12, e22) in additional_prd:
                    (D, weight, delta) = additional_prd[(e12.eid, e21.eid)][0]
                    if not self.graph.is_connected(first_rectangle.e2.v2, second_rectangle.e1, 10):
                        count_correct_scaffolds += 1
                    if len(first_rectangle.e2.v2.out) != 0 or len(second_rectangle.e1.v1.inn) != 0:
                        continue
                    used_paires.add((e12.eid, e21.eid))
                    count_incorrect_scaffolds += 1
                    if D - first_rectangle.e2.len > 0 and D - first_rectangle.e2.len < 100:
                        first_rectangle.e2.seq[-55:], "\n", second_rectangle.e1.seq[:55]
                        connect_edges.add((e1.eid, e2.eid))
                        max_eid = self.graph.max_eid
                        self.graph.add_edge(max_eid, e12.v2.vid, e21.v1.vid, self.graph.K + 3, max_eid + 1)
                        self.graph.add_edge(max_eid + 1, e21.conj.v2.vid, e12.conj.v1.vid, self.graph.K + 3, max_eid)
                        seq = first_rectangle.e2.seq[-self.graph.K:] + "NNN" + second_rectangle.e1.seq[:self.graph.K]
                        self.graph.add_seq(max_eid, seq)
                        self.graph.add_seq(max_eid + 1, utils.rc(seq))
                        path_1 = []
                        path_2 = []
                        used = set()
                        begin_path = False
                        start_offset = 0
                        for diag in e1.diagonals:
                            if e11 == diag.rectangle.e2:
                                begin_path = True
                            if begin_path and diag.rectangle.e2 not in used:
                                path_1.append(diag.rectangle.e2)
                                used.add(diag.rectangle.e2)
                        path_1.append(self.graph.es[max_eid])
                        if e1.diagonals[-1].rectangle.e2.len <= e1.diagonals[-1].offsetc:
                            path_1 = path_1[1:]
                            start_offset = 0
                        else:
                            start_offset = e1.diagonals[-1].offsetc
                        path_2.append(self.graph.es[max_eid])
                        path_2.append(e2.diagonals[0].rectangle.e1)
                        used = set()
                        for diag in e2.diagonals:
                            if e22 == diag.rectangle.e1:
                                break
                            if diag.rectangle.e1 not in used:
                                path_2.append(diag.rectangle.e1)
                                used.add(diag.rectangle.e1)

        self.test_utils.logger.info("count_correct_scaffolds " + str(count_correct_scaffolds) + " " + str(
            count_incorrect_scaffolds) + " " + str(len(used_paires)) + "\n")
        return connect_edges

    def add_rectangles_by_path(self, path1, path2, start_offset):
        path_len = 0
        for p in path1:
            path_len += p.len
            # path_len -= start_offset
        first_shift = start_offset
        second_shift = 0
        pos_first_path = 0
        pos_second_path = 0
        first_len = first_shift

        while first_len < path_len:
            ed1 = path1[pos_first_path]
            ed2 = path2[pos_second_path]
            rectangle = Rectangle(ed1, ed2)
            rectangle.add_diagonal(self.d, self.d + first_shift - second_shift)
            rect_diag = rectangle.get_closest_diagonal(self.d + first_shift - second_shift)
            self.add_diagonal_and_conj(rect_diag)
            if ed2.len - second_shift < ed1.len - first_shift:
                pos_second_path += 1
                first_shift += ed2.len - second_shift
                first_len += ed2.len - second_shift
                second_shift = 0
            elif ed1.len - first_shift < ed2.len - second_shift:
                pos_first_path += 1
                first_len += ed1.len - first_shift
                second_shift += ed1.len - first_shift
                first_shift = 0
            else:
                first_len += ed1.len - first_shift
                pos_second_path += 1
                pos_first_path += 1
                first_shift = 0
                second_shift = 0

    def edges_expand(self, L):
        should_connect = dict()
        for edge_id, edge in self.es.items():
            if edge.length() > L:
                self.edge_expand(edge, should_connect, L)
        to_delete = set()
        for edge_id, path in should_connect.items():
            if path[-1].eid in should_connect:
                to_delete.add(edge_id)
        for eid in to_delete:
            del should_connect[eid]
        conj_should_connect = dict()
        for edge_id, path in should_connect.items():
            conj_path = []
            i = len(path) - 1
            while i >= 0:
                conj_path.append(path[i].conj)
                i -= 1
            conj_should_connect[path[0].conj.eid] = conj_path
            conj_should_connect[path[0].eid] = path
        return conj_should_connect

    def edge_expand(self, begin_edge, should_connect, L):
        second_edges = []
        prev_first = None
        first = begin_edge.diagonals[0].rectangle.e2
        for diag in begin_edge.diagonals:
            if diag.rectangle.e1 == first:
                prev_first = first
                while prev_first.eid == first.eid:
                    if len(second_edges) == 0:
                        first == None
                        break
                    first = second_edges.pop(0)
            second_edges.append(diag.rectangle.e2)
        if len(second_edges) == 0:
            return
        v2 = begin_edge.v2
        paths = self.expand(v2, second_edges, first, prev_first, [begin_edge], [])
        if len(paths) == 1:
            should_add = False
            path = paths[0]
            for e in path:
                if e != begin_edge and e.length() > L:
                    should_add = True
            if should_add:
                should_connect[begin_edge.eid] = path
        return

    def expand(self, v2, second_edges, first, prev_first, path, paths):
        extend_edges = []
        end_edges = []
        for next_edge in v2.out:
            extend = self.can_expand(next_edge, second_edges, first, prev_first)
            if not extend:
                continue
            if len(extend[1]) == 0 and not extend[2]:
                new_path = list(path)
                new_path.append(next_edge)
                paths.append(new_path)
            elif len(extend[1]) >= 1:
                extend_edges.append(extend)
        for next_extend in extend_edges:
            (n_edge, n_second_edges, n_first, n_prev_first) = next_extend
            new_path = list(path)
            new_path.append(n_edge)
            self.expand(n_edge.v2, n_second_edges, n_first, n_prev_first, new_path, paths)
        return paths


    def can_expand(self, edge, second_edges, first, prev_first):
        second_edges = list(second_edges)
        for diag in edge.diagonals:
            if diag.rectangle.e1 != first and diag.rectangle.e1 != prev_first:
                return None
  
            if diag.rectangle.e1 == first:
                prev_first = first
                while prev_first.eid == first.eid:
                    if len(second_edges) == 0:
                        return (edge, second_edges, None, prev_first)
                    first = second_edges.pop(0)
        return (edge, second_edges, first, prev_first)


    def __try_delete_bv(self, v):
        if len(v.out) == 0 and len(v.inn) == 0 and v.key in self.vs:
            del self.vs[v.key]

    def check_tips(self, K ):
        tips = set()
        for bv in self.vs.itervalues():
            if len(bv.inn) == 1 and len(bv.out) == 0 and len(bv.inn[0].get_seq(K, self.d)) < 3 * self.d and bv.inn[
                                                                                                            0].v1.vid != bv.vid:
               
                tips.add(bv.inn[0])
        self.delete_tips(K, tips)

    def delete_tips(self, K, bes):
        for be in bes:
            self.__remove_bedge__(be)
            if (be.eid != be.conj.eid):
                self.__remove_bedge__(be.conj)
        self.condense()

    def path_is_subpath(self, path, sub_path):
        cur_edge_index = 0
        cur_edge_eid  = sub_path[cur_edge_index]
        for edge in path:
            if edge.eid == cur_edge_eid:
                cur_edge_index += 1
                if cur_edge_index == len(sub_path):
                    return True
            else:
                cur_edge_index = 0
            cur_edge_eid = sub_path[cur_edge_index]
        return False
            
    def choose_best_path(self, all_paths, long_eid1, long_eid2, begs):
        if len(all_paths) < 2:
            return all_paths[0]
        real_first_paths = []
        real_second_paths = []
        for be in begs:
            path_first = []
            path_second = []
            for diag in be.diagonals:
                e1_eid = diag.rectangle.e1.eid
                e2_eid = diag.rectangle.e2.eid
                if e1_eid == long_eid1:
                    path_first = []
                if e2_eid == long_eid1:
                    path_second = []
                if e1_eid == long_eid2:
                    path_first.append(e1_eid)
                    break
                if len(path_first) == 0 or path_first[-1] != e1_eid:
                    path_first.append(e1_eid)
                if len(path_second) == 0 or path_second[-1] != e2_eid:
                    path_second.append(e2_eid)
            if len(path_first) > 1:
                real_first_paths.append(path_first)
            if len(path_second) > 1:
                real_second_paths.append(path_second)
        best_count_real_paths = 0
        best_path = all_paths[0]
        for path in  all_paths:
            count_real_paths = 0
            path.append(self.graph.es[long_eid2])
            for real_path in real_first_paths:
                if self.path_is_subpath(path, real_path):
                    count_real_paths += 1
            for real_path in real_second_paths:
                if self.path_is_subpath(path, real_path):
                    count_real_paths += 1
            if count_real_paths > best_count_real_paths :
                best_count_real_paths = count_real_paths
                best_path = path
        return best_path[:-1]


    def delete_missing_loops(self, DG_loops, K, L, threshold):
        if DG_loops == None:
            return
        begs_related_to_loop = defaultdict(set)
        begin_loops = dict()
        end_loops = dict()
        for k, be in self.es.items():
            for diag in be.diagonals:
                for eeid1, (long_eid1, long_eid2, busheids, path, visited_vs, all_paths) in DG_loops.items():
                    rect = diag.rectangle
                    if not set([rect.e1.eid, rect.e2.eid]) <= busheids:
                        continue
                    
                    if rect.e1.eid == long_eid1 and rect.e2.eid == long_eid1:
                        begin_loops[long_eid1] = (diag, be)
                    
                    if rect.e1.eid == long_eid2 and  rect.e2.eid == long_eid2:
                        end_loops[long_eid1] = (diag, be)                        
                    
                    begs_related_to_loop[eeid1].add(be)

        diag_to_add = []
        for eid, begs in begs_related_to_loop.items():
            (long_eid1, long_eid2, busheids, path, visited_vs, all_paths) = DG_loops[eid]
            path = self.choose_best_path(all_paths, long_eid1, long_eid2, begs)
            if len(begs) < 2:
                continue
            if eid not in begin_loops or eid not in end_loops:
                continue
            begin_diag = begin_loops[eid][0]
            if begin_diag.rectangle.e1.len < 2*self.d or len(path) < 1 or path[0] != begin_diag.rectangle.e1 or path[0].len < 2 * self.d:
                self.logger.info("BAD CASE " + str( eid)+ " " + str( begin_diag.rectangle.e1.len) + " " +str((long_eid1, long_eid2, busheids, path, visited_vs))) 
                continue
            end_diag = end_loops[eid][0]
            path.append(end_loops[eid][0].rectangle.e1)
            first_shift = begin_diag.offseta
            second_shift = begin_diag.offsetb
            
            path_len = reduce(lambda s, e: s + e.len, path, 0)           
            rectangles = []
            diags = []
            pos_first_path = 0
            pos_second_path = 0
            first_len = first_shift
            while first_len < path_len and pos_second_path < len(path):
                ed1 = path[pos_first_path]
                ed2 = path[pos_second_path]
                rectangle = Rectangle(ed1, ed2)
                rectangle.add_diagonal(self.d, self.d + first_shift - second_shift)
                rect_diag = rectangle.get_closest_diagonal(self.d + first_shift - second_shift)
                rectangles.append(rectangle)
                diags.append(rect_diag)
                if ed2.len - second_shift < ed1.len - first_shift:
                    pos_second_path += 1
                    first_shift += ed2.len - second_shift
                    first_len += ed2.len - second_shift
                    second_shift = 0
                elif ed1.len - first_shift < ed2.len - second_shift:
                    pos_first_path += 1
                    first_len += ed1.len - first_shift
                    second_shift += ed1.len - first_shift
                    first_shift = 0
                else:
                    first_len += ed1.len - first_shift
                    pos_second_path += 1
                    pos_first_path += 1
                    first_shift = 0
                    second_shift = 0
            
            diag_to_add.extend(diags)
            
            for bedge in begs:
                self.__remove_bedge_and_conj(bedge)
             
                for diag in bedge.diagonals:
                    if diag.rectangle.e1.eid not in busheids or diag.rectangle.e2.eid not in busheids:
                        diag_to_add.append(diag)

                if begin_diag in bedge.diagonals:
                    for diag in bedge.diagonals:
                        if diag == begin_diag:
                            break
                        diag_to_add.append(diag)

                elif end_diag in bedge.diagonals:
                    bedge.diagonals.reverse()
                    for diag in bedge.diagonals:
                        if diag == end_diag:
                            break
                        diag_to_add.append(diag)

        for diag in diag_to_add:
            self.add_diagonal_and_conj(diag)


    #L - big contig's len > L
    #threshold - max deep in bfs
    def delete_loops(self, K, L, threshold):
        edges_to_delete = set()
        connected_paths = set()
        count_loop = 0
        rectangles_before_loop = []
        edges_before_loop = set()
        for k, be in self.es.items():
            if len(be.get_seq(K, self.d)) > L:
                found = self.__find_loops(be, K, threshold, L, edges_to_delete, connected_paths)
                if found:
                    count_loop += 1
                    rectangles_before_loop.append([diag.rectangle.e1.eid for diag in be.diagonals])
                    for diag in be.diagonals:
                        edges_before_loop.add(diag.rectangle.e1.eid)
                        edges_before_loop.add(diag.rectangle.e2.eid)
                        #for e in rectangles_before_loop[-1]:
                        #  edges_before_loop.add(e)
        for edge in edges_to_delete:
            if edge.eid not in connected_paths and edge.conj.eid not in connected_paths:
                self.__remove_bedge_and_conj(edge)
        return edges_before_loop

    def __find_loops(self, be, K, threshold, L, edges_to_delete, connected_path):
        result = self.find_all_loops(be, threshold, L, dict())
        if not result:
            return False
        long_end = self.es[result[1]]
        end_long_edge_id = long_end.v2.vid
        begin_long_edge_id = self.es[result[0]].v1.vid
        visited_vs = result[4]
        for bv in visited_vs:
            for e in bv.out:
                if e.v2.vid != end_long_edge_id:
                    edges_to_delete.add(e)
            for e in bv.inn:
                bv_begin = e.v1
                if bv_begin.vid != begin_long_edge_id:
                    edges_to_delete.add(e)
        path = self.get_paths(be.v2, long_end, None, threshold)[0]
        for e in path:
            connected_path.add(e.eid)
        return True

    def print_about_edges(self, eids, K):
        print "All about edge", eids
        for k, be in self.es.items():
            if be.eid in eids:
                info = str(be.eid) + " len " + str(len(be.get_seq(K, self.d))) + " "
                v1 = be.v1
                v2 = be.v2
                info += " v1 " + str(v1.vid) + " inn "
                for e in v1.inn:
                    info += str(e.eid) + " (" + str(e.v1.vid) + "," + str(e.v2.vid) + ") "
                info += " out "
                for e in v1.out:
                    info += str(e.eid) + " (" + str(e.v1.vid) + "," + str(e.v2.vid) + ") "
                info += " v2 " + str(v2.vid) + " inn "
                for e in v2.inn:
                    info += str(e.eid) + " (" + str(e.v1.vid) + "," + str(e.v2.vid) + ") "
                info += " out "
                for e in v2.out:
                    info += str(e.eid) + " (" + str(e.v1.vid) + "," + str(e.v2.vid) + ") "
                info += "\n"
                print info

    def save(self, outpath, K):
        eid = 0
        left = open(os.path.join(outpath, "rectangle_paired_info_1.fasta"), "w")
        right = open(os.path.join(outpath, "rectangle_paired_info_2.fasta"), "w")
        for k, be in self.es.items():
            (seq1, seq2) = be.get_paired_seq(K, self.d)
            seq2 = seq2.strip('N')
            seq1 = seq1.strip('N')
            left.write(">" + str(eid) + "/1\n" + seq1 + "\n")
            right.write(">" + str(eid) + "/2\n" + seq2 + "\n")
            eid += 1
        left.close()
        right.close()

    def __get_bvertex(self, key):
        if key in self.vs:
            #print key
            bv = self.vs[key]
            #bv.key = bv.key.join_with(key) # transitive closure
        else:
            self.vs[key] = BVertex(key)
        return self.vs[key]

    def __is_bvertex(self, key):
        return key in self.vs

    def __add_bedge(self, diag):
        v1 = self.__get_bvertex(diag.key1)
        v2 = self.__get_bvertex(diag.key2)
        be = BEdge(v1, v2, diag)
        self.es[be.eid] = be
        return be

    def add_diagonal(self, diag):
        if diag in self.diagonals:
            return
       
        be = self.__add_bedge(diag)
        conj = be
        if diag.conj != diag:
            conj = self.__add_bedge(diag.conj)
            #print len(self.es), len(self.vs)
        self.diagonals.add(diag)
        self.diagonals.add(diag.conj)
       # print diag.rectangle.e1.eid, diag.rectangle.e2.eid, diag.D, diag.conj.rectangle.e1.eid, diag.conj.rectangle.e2.eid, diag.conj.D
       # print diag.conj.rectangle.e1.eid, diag.conj.rectangle.e2.eid, diag.conj.D, diag.rectangle.e1.eid, diag.rectangle.e2.eid, diag.D
        conjugate(be, conj)
        conjugate(be.v1, conj.v2)
        conjugate(be.v2, conj.v1)
        return (be, conj)

    def add_diagonal_and_conj(self, diag):
        if diag in self.diagonals:
            return 
        rect = diag.rectangle
        rect_conj = Rectangle(rect.e2.conj, rect.e1.conj)
        conjugate(rect, rect_conj)
        D = diag.D - diag.rectangle.e1.len + diag.rectangle.e2.len
        if experimental.filter == experimental.Filter.pathsets:
            pathset = diag.pathset.conj()
        else:
            pathset = None
        rect_conj.add_diagonal(self.d, D, pathset)
        diag_conj = rect.conj.diagonals[D, pathset]
        conjugate(diag, diag_conj)

        return self.add_diagonal(diag)

    def __join_biedges(self, be1, be2):
        ## u ---be1---> v ---be2---> w
        ## z <--be4---- y <--be3---- x
        ## transforms to:
        ## u --------beA--------> w
        ## z <-------beB--------- x
        be3 = be2.conj
        be4 = be1.conj
        u, v, w = be1.v1, be1.v2, be2.v2
        x, y, z = be3.v1, be3.v2, be4.v2
        assert be1.v2 == be2.v1
        assert 1 == len(v.inn) == len(v.out) == len(y.inn) == len(y.out), (
        be1.eid, be2.eid, len(v.inn), len(v.out), len(y.inn), len(y.out))
        if be1 == be3 or be2 == be4:
            return
        assert be1 != be3, "=> (v == y) => (in-degree(v) > 1)"
        assert be2 != be4, "=> (v == y) => (out-degree(v) > 1)"

        if be1 == be4 and be2 == be3:
            assert z == v == x
            assert u == y == w
            #assert False
            return # TODO: think how to condense better, rare case!!!!

        if be2 == be3: # loop on the right: be1->be2=be3->be4
            assert v == x
            assert y == w
            beA = BEdge(u, z, None)
            beA.diagonals = be1.diagonals + be2.diagonals + be4.diagonals
            first_connect = self.test_utils.should_join(be1.diagonals[-1], be2.diagonals[0])
            second_connect = self.test_utils.should_join(be2.diagonals[-1], be4.diagonals[0])
            if first_connect:
                self.test_utils.join_correct += 1
            else:
                self.test_utils.join_incorrect += 1
            if second_connect:
                self.test_utils.join_correct += 1
            else:
                self.test_utils.join_incorrect += 1
            conjugate(beA, beA)
            self.es[beA.eid] = beA
            u.out.remove(be1)
            w.inn.remove(be2)
            z.inn.remove(be4)
            del self.es[be1.eid]
            del self.es[be2.eid]
            del self.es[be4.eid]
        elif be1 == be4: # loop on the left: be3->be1=be4->be2
            assert u == y
            assert z == v
            beA = BEdge(x, w, None)
            beA.diagonals = be3.diagonals + be1.diagonals + be2.diagonals
            first_connect = self.test_utils.should_join(be3.diagonals[-1], be1.diagonals[0])
            second_connect = self.test_utils.should_join(be1.diagonals[-1], be2.diagonals[0])
            if first_connect:
                self.test_utils.join_correct += 1
            else:
                self.test_utils.join_incorrect += 1
            if second_connect:
                self.test_utils.join_correct += 1
            else:
                self.test_utils.join_incorrect += 1
            conjugate(beA, beA)
            self.es[beA.eid] = beA
            u.out.remove(be1)
            w.inn.remove(be2)
            x.out.remove(be3)
            del self.es[be1.eid]
            del self.es[be2.eid]
            del self.es[be3.eid]
        else: # most usual case
            be_set = set()
            be_set.add(be1)
            be_set.add(be2)
            be_set.add(be3)
            be_set.add(be4)
            assert len(be_set) == 4, (be1, be2, be3, be4)  #all different
            six_set = set()
            six_set.add(u)
            six_set.add(v)
            six_set.add(w)
            six_set.add(x)
            six_set.add(y)
            six_set.add(z)
            """if u == w:
                assert z == x
                assert len(six_set) == 4, (u, v, w, x, y, z) # same ends, ok
            elif u == x:
                assert z == w
                assert len(six_set) == 4, (u, v, w, x, y, z) # conjugated ends, ok
            else:
                assert len(six_set) == 6, (u, v, w, x, y, z) # all different
                # TODO: check (x == u and w == z)"""
            beA = BEdge(u, w, None)
            beA.diagonals = be1.diagonals + be2.diagonals
            if self.test_utils:
                first_connect = self.test_utils.should_join(be1.diagonals[-1], be2.diagonals[0])
                second_connect = self.test_utils.should_join(be3.diagonals[-1], be4.diagonals[0])
                if first_connect:
                    self.test_utils.join_correct += 1
                else:
                    self.test_utils.join_incorrect += 1
                if second_connect:
                    self.test_utils.join_correct += 1
                else:
                    self.test_utils.join_incorrect += 1
            beB = BEdge(x, z, None)
            beB.diagonals = be3.diagonals + be4.diagonals
            conjugate(beA, beB)
            self.es[beA.eid] = beA
            self.es[beB.eid] = beB
            u.out.remove(be1)
            w.inn.remove(be2)
            x.out.remove(be3)
            z.inn.remove(be4)
            del self.es[be1.eid]
            del self.es[be2.eid]
            del self.es[be3.eid]
            del self.es[be4.eid]

        v.inn, v.out = [], []
        y.inn, y.out = [], []
        self.vs.pop(v.key)
        self.vs.pop(y.key)

    def condense(self):
        l = len(self.vs)
        for bv in self.vs.values(): # copy because can be: "Set changed size during iteration"
            if len(bv.inn) == 1 and len(bv.out) == 1 and (bv.inn[0] != bv.out[0]):
                self.__join_biedges(bv.inn[0], bv.out[0])
                # self.__check()
        self.logger.info("Condensed %d bi-vertices (out of %d). %d bi-vertices left." % (l - len(self.vs), l,
                                                                                         len(self.vs)))

    def project(self, outpath, is_sc, is_careful):
        g = graph.Graph()
        count = 0
        for key, be in self.es.items():
            #print count + 1 
            count += 1
            # v ---------be--------> w
            # y <-----be.conj------- x
            v, w = be.v1, be.v2
            x, y = be.conj.v1, be.conj.v2
            #assert v != w
            #assert v != x
            seq = be.get_seq_for_contig(self.graph.K, self.d, is_sc, is_careful)
            cvr = be.get_cvr()
            g.add_vertex(v.vid, y.vid)
            g.add_vertex(y.vid, v.vid)
            g.add_vertex(w.vid, x.vid)
            g.add_vertex(x.vid, w.vid)
            g.add_edge(be.eid, v.vid, w.vid, len(seq) - self.graph.K, be.conj.eid)
            g.add_seq(be.eid, seq)
            g.add_cvr(be.eid, cvr)
        g.update_K()
        maxv = BVertex.vid
        maxe = BEdge.eid
        taken = set()
        for diag in self.diagonals:
            taken.add(diag.rectangle.e1)
            taken.add(diag.rectangle.e2)
        for e in self.graph.es.itervalues():
            if e not in taken:
                # v ---e1---> w
                # x <--e2---- y
                assert e.conj not in taken
                e1 = e
                e2 = e.conj
                v = e1.v1.vid + maxv
                w = e1.v2.vid + maxv
                y = e2.v1.vid + maxv
                x = e2.v2.vid + maxv
                g.add_vertex(v, x)
                g.add_vertex(x, v)
                g.add_vertex(w, y)
                g.add_vertex(y, w)
                seq = e1.seq
                g.add_edge(e1.eid + maxe, v, w, len(seq) - self.graph.K, e2.eid + maxe)
                g.add_seq(e1.eid + maxe, seq)
        return g

    def __check(self):
        for eid, edge in self.es.items():
            for this, next in itertools.izip(edge.diagonals, edge.diagonals[1:]):
                assert this.key2 == next.key1, (this, "->", next)

    
