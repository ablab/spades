import logging
import sys 
import math
import os
import utils
import saveparser
from diagonal import Diagonal
from utils import conjugate
import experimental

def makelogger(logfilename):
  logger = logging.getLogger('debug')
  logger.setLevel(logging.DEBUG)
  fh = logging.FileHandler(logfilename, mode='w')
  fh.setLevel(logging.DEBUG)
  ch = logging.StreamHandler()
  ch.setLevel(logging.ERROR)
  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
  fh.setFormatter(formatter)
  ch.setFormatter(formatter)
  logger.addHandler(fh)
  logger.addHandler(ch)



def parse_ref_info(fin):
  ref_info = dict()
  for line in fin:
    block = line.split()
    edge_id = int(block[0].strip())
    ref_begin = int(block[1].strip())
    ref_end = int(block[2].strip())
    contig_begin = int(block[3].strip())
    contig_end = int(block[4].strip())
    if edge_id < 0:
      edge_id = -edge_id
      if edge_id not in ref_info:
        ref_info[edge_id] = []
      ref_info[edge_id].append(-(ref_begin - contig_begin))
    else:
      if contig_end > contig_begin:
        if edge_id not in ref_info:
          ref_info[edge_id] = []
        ref_info[edge_id].append(ref_begin - contig_begin)
  #print "ref_info", len(ref_info)
  return ref_info

class TestUtils(object):
  def __init__(self, edge_aligned_file, log_file_name):
    makelogger(log_file_name)
    self.logger = logging.getLogger('debug')
    self.has_ref_info = True  
    if edge_aligned_file:
      try:
        with open(edge_aligned_file) as f: 
           self.ref_info = parse_ref_info(open(edge_aligned_file))
      except IOError as e:
           self.ref_info = dict()
           self.has_ref_info = False
    else:
      self.ref_info = dict()
      self.has_ref_info = False
    self.UNALIGNED_DIAG = 0
    self.TRUE_DIAG = 1
    self.FALSE_DIAG = 2
    self.unaligned = 0
    self.not_true_diags = 0
    self.true_diags = 0
    self.join_correct = 0
    self.join_incorrect = 0
    self.join_unaligned = 0
    self.similar_diags = dict()
  
  def __get_ref_info(self, eid):
    if eid in self.ref_info:
      return self.ref_info[eid]
    return []

  def add_to_diags(self, diag):
    rect = diag.rectangle
    e1 = rect.e1
    e2 = rect.e2
    D = diag.D
    if (e1, e2) not in self.similar_diags:
      self.similar_diags[(e1, e2)] = []
    self.similar_diags[(e1, e2)].append(D)

  def print_similar_diags(self, coef):
    if not self.has_ref_info:
      return
    res = []
    for i in range(coef):
      res.append([])
    for (e1, e2), list_of_D in self.similar_diags.items():
      list_of_D.sort()
      for i in range(len(list_of_D)):
        for j in range(1, coef + 1):
          """while begin >= 0 and list_of_D[i] - list_of_D[begin] == j:
            res[j-1].append((e1,e2, list_of_D[i], list_of_D[begin]))
            begin = begin -1"""
          end = i+1
          while end < len(list_of_D) and list_of_D[end] - list_of_D[i] < j:
            end += 1
          while end < len(list_of_D) and list_of_D[end] - list_of_D[i] == j:
            res[j-1].append((e1,e2, list_of_D[i], list_of_D[end]))
            end += 1
    self.logger.info(res)
    for i in range(coef):
      self.logger.info(str(i + 1) + " " + str(len(res[i])))


  def is_true_diagonal(self, diag):
    rect = diag.rectangle
    e1 = rect.e1
    e2 = rect.e2
    D = diag.D
    if e1.eid == e2.eid:
      return self.TRUE_DIAG
    e1_ref_infos = self.__get_ref_info(e1.eid)
    e2_ref_infos = self.__get_ref_info(e2.eid)
    if len(e1_ref_infos) == 0 or len(e2_ref_infos) == 0:
      self.unaligned += 1
      #return False
      return self.UNALIGNED_DIAG
    true_diag = False
    for e1_ref_info in e1_ref_infos:
      for e2_ref_info in e2_ref_infos:
        dist = abs(e2_ref_info - e1_ref_info)
        if dist == D or dist == D + 1 or dist == D - 1 or dist == D + 2 or dist == D - 2:
          true_diag = True
          break
      if true_diag:
        break
    if true_diag:
      self.true_diags += 1
    else:
      self.not_true_diags += 1
    ###print "true_diag", true_diag
    """if not true_diag:
      print "edges", e1.eid, e2.eid
      print "e1",  e1_ref_infos
      print "e2", e2_ref_infos
      print "D", D

    print "diag support", diag.support(), "len", diag.offsetc- diag.offseta"""
    
    return self.TRUE_DIAG if true_diag else self.FALSE_DIAG

  def should_join(self, diag1, diag2):
     rect1 = diag1.rectangle
     rect2 = diag2.rectangle
     e1_1 = rect1.e1
     e1_2 = rect1.e2
     e2_1 = rect2.e1
     e2_2 = rect2.e2
     offseta_2 = diag2.offseta
     offsetb_2 = diag2.offsetb
     offsetc_1 = diag1.offsetc
     offsetd_1 = diag1.offsetd
     if len(self.__get_ref_info(e1_1.eid)) == 0 or len(self.__get_ref_info(e1_2.eid)) == 0 or len(self.__get_ref_info(e2_1.eid) ) == 0 or len(self.__get_ref_info(e2_2.eid)) == 0:
      self.join_unaligned += 1
      return True

     for e11_ref_info in self.__get_ref_info(e1_1.eid):
      for e12_ref_info in self.__get_ref_info(e1_2.eid):
        for e21_ref_info in self.__get_ref_info(e2_1.eid):
          for e22_ref_info in self.__get_ref_info(e2_2.eid):
            if e11_ref_info*e21_ref_info < 0 or e12_ref_info * e22_ref_info < 0:
             # print "join different strand"
              continue
            #print "join",  e11_ref_info - e21_ref_info, abs(e11_ref_info) + offsetc_1 -(  abs(e21_ref_info) + offseta_2), e12_ref_info - e22_ref_info, abs(e12_ref_info) + offsetd_1 -(abs(e22_ref_info) + offsetb_2)
            if abs(e11_ref_info) + offsetc_1 == abs(e21_ref_info) + offseta_2 and abs(e12_ref_info) + offsetd_1 == abs(e22_ref_info) + offsetb_2:
              return True
     #print offsetc_1, offseta_2, offsetd_1, offsetb_2
     #print self.__get_ref_info(e1_1.eid)
     #print self.__get_ref_info(e1_2.eid)
     #print self.__get_ref_info(e2_1.eid)
     #print self.__get_ref_info(e2_2.eid)
     return False
  
  def stat(self):
    self.logger.info("unaligned " + str(self.unaligned) +  " true_diags " + str(self.true_diags) +  " not_true_diags " + str(self.not_true_diags) +  " join correct " + str(self.join_correct) +  " join incorrect"  + str(self.join_incorrect) + " join unaligned " + str(self.join_unaligned))

  def have_ref(self, e_first, e_prev, e_cur, pos):
    if len(self.__get_ref_info(e_first.eid)) == 0 or len(self.__get_ref_info(e_prev.eid)) == 0 or len(self.__get_ref_info(e_cur.eid)) == 0 :
      return False
    for e_first_ref in self.__get_ref_info(e_first.eid):
      for e_prev_ref in self.__get_ref_info(e_prev.eid):
        for e_cur_ref in self.__get_ref_info(e_cur.eid):
          if e_first_ref * e_prev_ref < 0 or e_prev_ref * e_cur_ref < 0:
            continue
            if abs(e_cur_ref) == abs(e_prev_ref) + e_prev.len and abs(e_cur_ref) == abs(e_first_ref) + pos + e_first.len:
	      return True
    return False

  def count_of_genomic_pathes(self, e1, e2, D):
    if D == 0 or len(self.__get_ref_info(e1.eid)) == 0 or len(self.__get_ref_info(e2.eid)) == 0:
      return 0
    vertex_dist = int(D - e1.len)
    limit = vertex_dist + 1			
    ls = [{} for _ in xrange(limit)]
    ls[0][e1] = 1
    for pos in xrange(limit):
      for e, cnt in ls[pos].iteritems():
	v = e.v2
        for e_out in v.out:
	  pos2 = pos + e_out.len
 	  if pos2 < limit and self.have_ref(e1, e, e_out, pos, ptc):
            ls[pos2][e_out] = ls[pos2].get(e_out, 0) + 1
											
      count = 0
      for e, cnt in ls[vertex_dist].iteritems():
	if e.v2 == e2.v1 and self.have_ref(e1, e, e2, vertex_dist, ptc):
	  count += 1
    return count
