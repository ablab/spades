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
import utils

DEEP_THRESHOLD = 10

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
    length = 0
    for diag in self.diagonals:
      length += diag.offsetc - diag.offseta
    return length
  
  def get_begin_seq(self, K, d, is_sc):
    (seq1, seq2) = self.get_paired_seq(K, d)
    seq = ''.join(map(lambda x, y: x if x != 'N' else (y if y != 'N' else 'N'), seq1, seq2)).strip('N')
    seq = seq.split(self.get_midle_seq())[0]
    return seq
   
  def get_end_seq(self, K, d, is_sc):
    (seq1, seq2) = self.get_paired_seq(K, d)
    seq = ''.join(map(lambda x, y: x if x != 'N' else (y if y != 'N' else 'N'), seq1, seq2)).strip('N')
    seq = seq.split(self.get_midle_seq())[1]
    return seq
   
  def get_midle_seq(self):
    seq = ""
    for diag in self.diagonals:
      seq += diag.rectangle.e1.seq[diag.offseta:diag.offsetc]
    return seq

  def get_seq_for_contig(self, K, d, is_sc):
    if is_sc:
      CUT_THRESHOLD = 2.0 #TODO: should take from histogram  
      CUT_LENGTH_THRESHOLD = 5.0
      MIN_LENGTH = 4 * d
    else:
      CUT_LENGTH_THRESHOLD = 5
      CUT_THRESHOLD = 0.0
      MIN_LENGTH = 0
    (seq1, seq2) = self.get_paired_seq(K, d)
    seq = ''.join(map(lambda x, y: x if x != 'N' else (y if y != 'N' else 'N'), seq1, seq2)).strip('N')
    first = self.diagonals[0]
    last = self.diagonals[-1]
    if len(seq1) > MIN_LENGTH:
      cur_len = 0
      diag_index = 0
      diag = self.diagonals[diag_index]
      can_add_begin = True
      while cur_len < d:
        if diag.offsetc - diag.offseta < CUT_LENGTH_THRESHOLD or (diag.support() < CUT_THRESHOLD):
          can_add_begin = False
          break
        diag_index += 1
        if diag_index == len(self.diagonals):
          cur_len = d
          continue
        diag = self.diagonals[diag_index]
        cur_len += diag.offsetc - diag.offseta
      
      cur_len = 0
      diag_index = len(self.diagonals) -1
      diag = self.diagonals[diag_index]
      can_add_end = True
      while cur_len < d:
        if diag.offsetc - diag.offseta < CUT_LENGTH_THRESHOLD or (diag.support() < CUT_THRESHOLD):
          can_add_end = False
          break
        diag_index -= 1
        if diag_index == -1:
          cur_len = d
          continue
        diag = self.diagonals[diag_index]
        cur_len += diag.offsetc - diag.offseta
      if can_add_end and can_add_begin:
        return first.rectangle.e1.seq[:first.offseta] + seq + last.rectangle.e2.seq[last.offsetd + K:] 
      if can_add_end:
        return seq + last.rectangle.e2.seq[last.offsetd + K:] 
      if can_add_begin:
        return first.rectangle.e1.seq[:first.offseta] + seq
    seq1 = cStringIO.StringIO()
    for this in self.diagonals:
      seq1.write(this.rectangle.e1.seq[this.offseta : this.offsetc])
    last = self.diagonals[-1]
    seq1.write(last.rectangle.e1.seq[last.offsetc : ])#last.offsetc + d])
    first = self.diagonals[0]
    seq1 = first.rectangle.e2.seq[:first.offsetb] + seq1.getvalue()[d:]#[first.offsetb - d:first.offsetb] + seq1.getvalue()[d:]  
    return seq1
  
  def get_seq(self,K,d):    
    (seq1, seq2) = self.get_paired_seq(K, d)
    seq = utils.seq_join(seq1, seq2).strip('N')
    return seq
  
  def get_paired_seq(self, K, d):
    seq1 = cStringIO.StringIO()
    seq2 = cStringIO.StringIO()
    seq2.write('N' * d)
    for this in self.diagonals:
      seq1.write(this.rectangle.e1.seq[this.offseta : this.offsetc])
      seq2.write(this.rectangle.e2.seq[this.offsetb : this.offsetd])
    last = self.diagonals[-1]
    seq1.write(last.rectangle.e1.seq[last.offsetc : last.offsetc + K])
    seq2.write(last.rectangle.e2.seq[last.offsetd : last.offsetd + K])
    seq1.write('N' * (len(seq2.getvalue())-len(seq1.getvalue())))
    return (seq1.getvalue(), seq2.getvalue())

  def get_cvr(self):
    cvr = 0.0
    sumlen = 0
    for this in self.diagonals:
      thiscvr = (this.rectangle.e1.cvr + this.rectangle.e2.cvr) * 0.5
      l =  this.offsetc - this.offseta
      cvr += thiscvr * l
      sumlen += l
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
        if (e12.eid, e21.eid) in additional_prd: #or (e11, e22) in additional_prd or (e11, e21) in additional_prd or (e12, e22) in additional_prd:
          (D, weight, delta) = additional_prd[(e12.eid,e21.eid)][0]
          if not self.graph.is_connected(first_rectangle.e2.v2, second_rectangle.e1, 10):
            count_correct_scaffolds +=1
          if len(first_rectangle.e2.v2.out) != 0 or len(second_rectangle.e1.v1.inn) != 0:
            continue
          used_paires.add((e12.eid, e21.eid))
          count_incorrect_scaffolds +=1
          if D - first_rectangle.e2.len > 0 and D - first_rectangle.e2.len < 100:
            print "SHOULD CONNECT", (e1.eid, e2.eid), (e12.eid, e21.eid), D - first_rectangle.e2.len, "\n", first_rectangle.e2.seq[-55:], "\n", second_rectangle.e1.seq[:55]
            connect_edges.add((e1.eid, e2.eid))
            max_eid = self.graph.max_eid
            self.graph.add_edge(max_eid, e12.v2.vid, e21.v1.vid, self.graph.K + 3, max_eid + 1)
            self.graph.add_edge(max_eid + 1, e21.conj.v2.vid, e12.conj.v1.vid, self.graph.K + 3, max_eid)
            seq = first_rectangle.e2.seq[-self.graph.K:] + "NNN" +  second_rectangle.e1.seq[:self.graph.K]
            self.graph.add_seq(max_eid, seq)
            self.graph.add_seq(max_eid + 1, utils.rc(seq))
            seq2 = second_rectangle.conj.e2.seq[-self.graph.K:] + "NNN" + first_rectangle.conj.e1.seq[:self.graph.K]
            assert seq2 == utils.rc(seq),"\n" +  seq2 + "\n" + utils.rc(seq)
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
                path_2.append(diag.rectange.e1)
                used.add(diag.rectangle.e1)
            print "path1", [e.eid for e in path_1] , "path2", [e.eid for e in path_2]  
            #self.add_rectangles_by_path(path_1, path_2, start_offset)
           
    self.test_utils.logger.info("count_correct_scaffolds " + str(count_correct_scaffolds) + " " + str(count_incorrect_scaffolds) + " " + str(len(used_paires)) + "\n")
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
          rectangle = Rectangle(ed1,ed2)
          rectangle.add_diagonal(self.d, self.d + first_shift - second_shift)
          rect_diag = rectangle.get_closest_diagonal(self.d + first_shift - second_shift) 
          self.add_diagonal_and_conj(rect_diag)
          print "ADD DIAGS", rect_diag
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
      i = len(path) -1
      while i >= 0:
        conj_path.append(path[i].conj)
        i -= 1
      conj_should_connect[path[0].conj.eid] = conj_path
      conj_should_connect[path[0].eid] = path
    for edge_id, path in conj_should_connect.items():
          print "NEW CONNECT PATHS",edge_id, [e.eid for e in path]
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
        else:
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
    v1s = set()
    v2s = set()
    tips = set()
    for bv in self.vs.itervalues():
        if len(bv.inn) == 1 and len(bv.out) == 0 and len(bv.inn[0].get_seq(K, self.d)) < 3 * self.d and bv.inn[0].v1.vid != bv.vid:
          edge =  bv.inn[0]
          if len(edge.diagonals) == 1:
            rect = edge.diagonals[0].rectangle
          v1s.add(bv)
          supp = 0
          for diag in edge.diagonals:
            supp += diag.support()
          tips.add(bv.inn[0])
    self.delete_tips(K,tips)

  def delete_tips(self, K, bes):
    for be in bes:
      self.__remove_bedge__(be)
      if (be.eid != be.conj.eid):
        self.__remove_bedge__(be.conj)
    self.condense()
  
  def delete_missing_loops(self, DG_loops, K, L, threshold):
    begs_related_to_loop = dict()
    begin_loops = dict()
    end_loops = dict()
    for eeid1, (long_eid1, long_eid2, busheids, path, visited_vs) in DG_loops.items():
      for k, be in self.es.items():
        for diag in be.diagonals:
          rect = diag.rectangle
          eids = [rect.e1.eid, rect.e2.eid ]
          if rect.e1.eid not in busheids or rect.e2.eid not in busheids:
            continue
          for eid in eids:
            if eid not in busheids:
              continue
            if rect.e1.eid == long_eid1:
              if rect.e2.eid == long_eid1: 
                begin_loops[long_eid1] = (diag, be)
            if rect.e2.eid == long_eid2:
              if rect.e1.eid == long_eid2:  
                end_loops[long_eid1] = (diag, be)
            if eeid1 not in begs_related_to_loop:
              begs_related_to_loop[eeid1] = set()
            begs_related_to_loop[eeid1].add(be)
    diag_to_add = set()
    for eid, begs in begs_related_to_loop.items():
      (long_eid1, long_eid2, busheids, path, visited_vs)  =  DG_loops[eid]
      if len(begs) < 2:
        continue
      if eid not in begin_loops or eid not in end_loops:
        print "not find begin_end"
        continue
      begin_diag = begin_loops[eid][0]
      end_diag = end_loops[eid][0]
      path.append(end_loops[eid][0].rectangle.e1)
      first_shift = begin_diag.offseta
      second_shift = begin_diag.offsetb
      path_len = 0
      for e in path:
        path_len += e.len
      rectangles = []
      diags = []
      pos_first_path = 0
      pos_second_path = 0
      first_len = first_shift
      while first_len < path_len and pos_second_path < len(path):
        ed1 = path[pos_first_path]
        ed2 = path[pos_second_path]
        rectangle = Rectangle(ed1,ed2)
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
      for diag in diags:
        diag_to_add.add(diag)
        
      for bedge in begs:
        self.__remove_bedge__(bedge)
        self.__remove_bedge__(bedge.conj)
        
        for diag in bedge.diagonals:
          if diag.rectangle.e1.eid not in busheids or diag.rectangle.e2.eid not in busheids:
            diag_to_add.add(diag)
        
        if begin_diag in bedge.diagonals:
          for diag in bedge.diagonals:
            if diag == begin_diag:
              break
            diag_to_add.add(diag)
          
        elif end_diag in bedge.diagonals:
          bedge.diagonals.reverse()
          for diag in bedge.diagonals:
            if diag == end_diag:
              break
            diag_to_add.add(diag)

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
        self.__remove_bedge__(edge)
        self.__remove_bedge__(edge.conj)
    return edges_before_loop

  def __find_loops(self, be, K, threshold, L, edges_to_delete, connected_path):
    result = self.find_all_loops(be, threshold, L)
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
    path = self.get_paths(be.v2, long_end, threshold)[0]
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
          info += str(e.eid) + " (" +str(e.v1.vid) + "," +  str(e.v2.vid) + ") "
        info += " out "
        for e in v1.out:
          info += str(e.eid) + " (" +str(e.v1.vid) + "," +  str(e.v2.vid) + ") "
        info += " v2 "  + str (v2.vid) + " inn "
        for e in v2.inn:
          info += str(e.eid) + " (" +str(e.v1.vid) + "," +  str(e.v2.vid) + ") "
        info += " out "
        for e in v2.out:
          info += str(e.eid) + " (" +str(e.v1.vid) + "," +  str(e.v2.vid) + ") "
        info += "\n" 
        print info

  def save(self, outpath, K):
    eid = 0
    left = open(os.path.join(outpath,"rectangle_paired_info_1.fasta"), "w")
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
      bv = self.vs.pop(key)
      bv.key = bv.key.join_with(key) # transitive closure
    else:
      bv = BVertex(key)
    self.vs[bv.key] = bv
    return bv

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
    conj = self.__add_bedge(diag.conj) if diag.conj != diag else be
    self.diagonals.add(diag)
    self.diagonals.add(diag.conj)
    conjugate(be, conj)
    conjugate(be.v1, conj.v2)
    conjugate(be.v2, conj.v1)
    return (be, conj)
    
  def add_diagonal_and_conj(self, diag):
    for old_diag in self.diagonals:
      if diag.rectangle.e1 == old_diag.rectangle.e1 and diag.rectangle.e2 == old_diag.rectangle.e2:
        if diag.D == old_diag.D:
          return
    rect = diag.rectangle
    rect_conj = Rectangle(rect.e2.conj, rect.e1.conj)
    conjugate(rect, rect_conj)
    D = diag.D - diag.rectangle.e1.len + diag.rectangle.e2.len       
    pathset = diag.pathset.conj() if experimental.filter == experimental.Filter.pathsets else None
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
        assert 1 == len(v.inn) == len(v.out) == len(y.inn) == len(y.out), (be1.eid, be2.eid, len(v.inn), len(v.out), len(y.inn), len(y.out))
        assert be1 != be3, "=> (v == y) => (in-degree(v) > 1)"
        assert be2 != be4, "=> (v == y) => (out-degree(v) > 1)"

        if be1 == be4 and be2 == be3:
            assert z == v == x
            assert u == y == w
            assert False
            return # TODO: think how to condense better, rare case

        if be2 == be3: # loop on the right: be1->be2=be3->be4
            assert v == x
            assert y == w
            beA = BEdge(u, z, None)
            beA.diagonals = be1.diagonals + be2.diagonals + be4.diagonals
            first_connect =  self.test_utils.should_join(be1.diagonals[-1], be2.diagonals[0])
            second_connect =  self.test_utils.should_join(be2.diagonals[-1],be4.diagonals[0]) 
            if first_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect += 1
            if second_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect +=1
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
            first_connect =  self.test_utils.should_join(be3.diagonals[-1], be1.diagonals[0])
            second_connect =  self.test_utils.should_join(be1.diagonals[-1],be2.diagonals[0]) 
            if first_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect += 1
            if second_connect:
              self.test_utils.join_correct += 1
            else:
              self.test_utils.join_incorrect +=1
            conjugate(beA, beA)
            self.es[beA.eid] = beA
            u.out.remove(be1)
            w.inn.remove(be2)
            x.out.remove(be3)
            del self.es[be1.eid]
            del self.es[be2.eid]
            del self.es[be3.eid]
        else: # most usual case
            assert len({be1, be2, be3, be4}) == 4, (be1, be2, be3, be4) # all different
            if u == w:
                assert z == x
                assert len({u, v, w, x, y, z}) == 4, (u, v, w, x, y, z) # same ends, ok
            elif u == x:
                assert z == w
                assert len({u, v, w, x, y, z}) == 4, (u, v, w, x, y, z) # conjugated ends, ok
            else:
                assert len({u, v, w, x, y, z}) == 6, (u, v, w, x, y, z) # all different
            # TODO: check (x == u and w == z)
            beA = BEdge(u, w, None)
            beA.diagonals = be1.diagonals + be2.diagonals
            if self.test_utils:
              first_connect =  self.test_utils.should_join(be1.diagonals[-1], be2.diagonals[0])
              second_connect =  self.test_utils.should_join(be3.diagonals[-1],be4.diagonals[0]) 
              if first_connect:
                self.test_utils.join_correct += 1
              else:
                self.test_utils.join_incorrect += 1
              if second_connect:
                self.test_utils.join_correct += 1
              else:
                self.test_utils.join_incorrect +=1
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

  def project(self, outpath, is_sc):
        log = open(os.path.join(outpath,"mis_log.txt"),"w")    
        g = graph.Graph()
        count_correct_rectangles = 0
        count_part_correct_rectangles = 0
        count_not_correct_rectanlges = 0
        count_part_unaligned_correct = 0
        count_unaligned = 0
        for key, be in self.es.items():
            # v ---------be--------> w
            # y <-----be.conj------- x
            v, w = be.v1, be.v2
            x, y = be.conj.v1, be.conj.v2
            correct_diag = 0
            unalign_diag = 0
            false_diag = 0
            for diag in be.diagonals:
              if self.test_utils:
                is_true_diag = self.test_utils.is_true_diagonal(diag)
                if is_true_diag == self.test_utils.TRUE_DIAG:
                  correct_diag +=1
                elif is_true_diag == self.test_utils.UNALIGNED_DIAG:
                  unalign_diag +=1
                else:
                  false_diag += 1
            if correct_diag > 0 and false_diag == 0 and unalign_diag == 0:
              count_correct_rectangles += 1
            elif correct_diag > 0 and false_diag == 0:
              count_part_unaligned_correct += 1
            elif correct_diag > 0:
              count_part_correct_rectangles += 1
            elif false_diag == 0:
              count_unaligned += 1
            elif false_diag > 0:
              count_not_correct_rectanlges +=1
            #assert be != be.conj
            #assert v != w
            #assert v != x
            seq = be.get_seq_for_contig(self.graph.K, self.d, is_sc)
            cvr = be.get_cvr()
            g.add_vertex(v.vid, y.vid)
            g.add_vertex(y.vid, v.vid)
            g.add_vertex(w.vid, x.vid)
            g.add_vertex(x.vid, w.vid)
            g.add_edge(be.eid, v.vid, w.vid, len(seq) - self.graph.K, be.conj.eid)
            g.add_seq(be.eid, seq)
            g.add_cvr(be.eid, cvr)
             
            log.write("\nmisassemble " + str(be.eid) + " "+ str(be.conj.eid)+ " "+ str(len(seq)))
            accum = 0
            for diag in be.diagonals:
                accum += diag.offsetc - diag.offseta
                log.write("\n" +  str(diag.offsetc - diag.offseta) + " " + str( accum) + " "+str(  diag.support()) + " diag.e1.len " +  str(diag.rectangle.e1.len) + " diag.e2.len " + str(diag.rectangle.e2.len)+ " e1.eid " + str(diag.rectangle.e1.eid) + " e2.eid " + str(diag.rectangle.e2.eid) )
        log.close()    
        if self.test_utils:
          self.test_utils.logger.info("count_correct_rectangles  = " + str(count_correct_rectangles) + "\ncount_part_unaligned_correct = " + str(count_part_unaligned_correct) + "\ncount_part_correct_rectangles  = " + str(count_part_correct_rectangles) + "\ncount_unaligned = " + str(count_unaligned) + "\ncount_not_correct = " + str(count_not_correct_rectanlges) + "\n\n") 
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

  def build_missing_rectangles(self, K, rectangles):
    return
    threshold = self.d
    self.test_utils.logger.info( "treshold " + str( threshold))
    count_ovelaps = 0
    count_miss_rect = 0
    count_miss_path = 0
    true_miss_path = 0
    count_overlaps = 0
    v1s = set()
    v2s = set()
    for bv in self.vs.itervalues():
        if len(bv.inn) == 1 and len(bv.out) == 0 and len(bv.inn[0].get_seq(K, self.d)) > 3 * self.d:
          v1s.add(bv)
        if len(bv.inn) == 0 and len(bv.out) == 1 and len(bv.out[0].get_seq(K, self.d)) > 3 * self.d:
          v2s.add(bv)
    assert len(v1s) == len(v2s) # because of rev-compl
    self.test_utils.logger.info("v1s.len "+ str( len(v1s)))
         
    all_paired_paths = []
    for v1 in v1s:
      be1 = v1.inn[0]
      diag1 = be1.diagonals[-1]
      for v2 in v2s:
        be2 = v2.out[0]
        if (be1.eid == be2.eid):
          continue
        diag2 = be2.diagonals[0]
        paths1 = abstract_graph.find_paths(diag1.rectangle.e1.v1, diag2.rectangle.e1.v1, diag1.rectangle.e1, threshold + diag1.rectangle.e1.len, DEEP_THRESHOLD)
        paths2 = abstract_graph.find_paths(diag1.rectangle.e2.v1, diag2.rectangle.e2.v1, diag1.rectangle.e2, threshold + diag1.rectangle.e2.len, DEEP_THRESHOLD)
        paired_paths = find_pair_paths(paths1, paths2, diag1, diag2)
        if len(paired_paths) != 0:
          all_paired_paths.append((paired_paths, diag1, diag2))
    self.test_utils.logger.info("all_paired_paths " + str( len(all_paired_paths)))
    can_find_one_path_more = True
    added_paths = []
    while can_find_one_path_more:
      the_best_path = None
      the_best_support = 0
      can_find_one_path_more = False

      for paired_paths in all_paired_paths:
        (best_support, best_len, best_rectangles, best_diags, best_path) = self.choose_best_path(paired_paths[0], rectangles, paired_paths[1], paired_paths[2], self.d, added_paths)
        if best_support/best_len >= 0.0 and best_support/best_len > the_best_support:
          the_best_support = best_support/best_len
          the_best_path = (best_support, best_len, best_rectangles, best_diags, best_path)
      if the_best_path:
        added_paths.append(the_best_path[-1])
        (best_support, best_len, best_rectangles, best_diags, best_path) = the_best_path
        can_find_one_path_more = True
        prev_diag = best_diags[0]
        true_path = True 
        for diag in best_diags[1:]:
          if prev_diag:
            should_connect = self.test_utils.should_join(prev_diag, diag)
            if not should_connect:
              true_path = False
            self.add_diagonal_and_conj(diag)
            is_true = self.test_utils.is_true_diagonal(diag)
            if not is_true:
              true_path = False
            count_miss_rect += 1
            prev_diag = diag
        count_miss_path += 1
        if true_path:
          true_miss_path += 1
    self.test_utils.logger.info( "count_overlap " + str( count_ovelaps) +  " count_miss_rect " + str( count_miss_rect) +  " count miss path " + str(count_miss_path) +  " true miss path " + str(true_miss_path))

  
    
  def choose_best_path(self, paired_paths, rectangeles_set, diag1, diag2, d, added_paths):
      best_support = 0
      best_len = 10000
      best_rectangles = []
      best_diags = []
      best_path = paired_paths[0]
      best_not_supported = 0
        
      for paired_path in paired_paths:
        (path1, path2, path_len) = paired_path
        if paired_path in added_paths:
          continue
        first_shift = diag1.offseta
        second_shift = diag1.offsetb
        path1.append(diag2.rectangle.e1)
        path2.append(diag2.rectangle.e2)
        rectangles = []
        diags = []
        not_supported = []
        path_support = 0
        pos_first_path = 0
        pos_second_path = 0
        first_len = first_shift
        make_less_N50 = False
        while not make_less_N50 and first_len < path_len + diag2.offseta:
          ed1 = path1[pos_first_path]
          ed2 = path2[pos_second_path]
          rectangle = Rectangle(ed1,ed2)
          rectangle.add_diagonal(d, d + first_shift - second_shift)
          rect_diag = rectangle.get_closest_diagonal(d + first_shift - second_shift) 
          if (not (rect_diag.key1 == diag1.key1 and rect_diag.key2 == diag1.key2) and not(rect_diag.key1 == diag2.key1 and rect_diag.key2 == diag2.key2)):
            can_use = [diag1.key1, diag1.key2, diag2.key1, diag2.key2]
            if (rect_diag.key1 in self.vs and rect_diag.key1 not in can_use) or  (rect_diag.key2 in self.vs and rect_diag.key2 not in can_use):
              make_less_N50 = True
              continue
          diags.append(rect_diag)
          rectangles.append(rectangle)
          rectangeles_set.use_prd_diag(rect_diag)
          #if rect_diag.prd_support < 0.00001 and (ed2.len > 10 and ed1.len > 10):
          #  make_less_N50 = True
          #  continue
          path_support += rect_diag.prd_support 
          if ed2.len - second_shift < ed1.len - first_shift:
            pos_second_path += 1
            first_shift += ed2.len - second_shift
            first_len += ed2.len - second_shift
            if rect_diag.prd_support < 0.000000001:
              not_supported.append(ed2.len - second_shift)
            second_shift = 0
          elif ed1.len - first_shift < ed2.len - second_shift:
            pos_first_path += 1
            first_len += ed1.len - first_shift
            second_shift += ed1.len - first_shift
            if rect_diag.prd_support < 0.000000001:
              not_supported.append(ed1.len - first_shift)
            first_shift = 0
          else:
            first_len += ed1.len - first_shift
            pos_second_path += 1
            pos_first_path += 1
            first_shift = 0
            second_shift = 0
        if not make_less_N50 and path_len > 1 and  path_support / path_len < 1000 and  path_support / path_len > best_support:
          best_support = path_support 
          best_len = path_len
          best_rectangles = rectangles
          best_diags = diags
          best_path = (path1, path2, path_len)
          best_not_supported = not_supported
      return (best_support,best_len, best_rectangles, best_diags, best_path)


def find_pair_paths(paths1, paths2, diag1, diag2):
  paired_paths = []
  for (path1, len1) in paths1:
    for (path2, len2) in paths2:
      if path1[0] != diag1.rectangle.e1 or path2[0] != diag1.rectangle.e2:
        continue
      if len1 - diag1.offseta + diag2.offseta == len2 - diag1.offsetb + diag2.offsetb: 
        paired_paths.append((path1, path2, len1))
  return paired_paths
  
  

