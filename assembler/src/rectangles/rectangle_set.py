import logging
from diagonal import Diagonal
from bigraph import BGraph
from utils import conjugate
from rectangle import Rectangle
from test_util import TestUtils
import pathsets
import experimental
import utils
import saveparser
import graph

class RectangleSet(object):

    def __init__(self, graph, d, test_utils = None, prd_file_name = None, first_prd_file_name = None, config = None):
        self.graph = graph
        self.d = d
        self.prd = dict()
        self.config = config
        self.additional_prd = dict() 
        if prd_file_name:
          self.__get_prd(prd_file_name)
        if first_prd_file_name:
          self.__get_additional_prd(first_prd_file_name)
        self.rectangles = {} # (e1, e2) -> Rectangle
        self.logger = logging.getLogger('rectangles')
        self.test_utils = test_utils
        self.not_used_prd_support = dict()
        
    def __get_prd(self, prd_file_name):
      for e1id, e2id, D, weight, delta in saveparser.prd(prd_file_name):
        self.prd[(e1id, e2id)] = (D, weight, delta + 2)
    
    def __get_additional_prd(self, prd_file_name):
      for e1id, e2id, D, weight, delta in saveparser.prd(prd_file_name):
        if (e1id, e2id) not in self.prd:
          e1= self.graph.es[e1id]
          e2 = self.graph.es[e2id]
          if len(e1.v2.out) != 0 or len(e2.v1.inn) != 0:
            continue
          e1_len = e1.len
          if D - e1_len > 0 and D - e1_len < 100:
            if (e1id, e2id) not in self.additional_prd:
              self.additional_prd[(e1id, e2id)] = []   
            self.additional_prd[(e1id, e2id)].append((D, weight, delta))


    def filter_without_prd(self):
      self.__build_from_graph()
      self.__conjugate()

    def filter(self, prd_filename, config):
        self.__build_from_graph()
        self.__conjugate()
        self.__use_prd(prd_filename, config)
        self.__rank()

    def pathsets(self, pst_filename):
        self.__build_from_pst(pst_filename)
        self.__conjugate()
        self.__all()

    def percentiles(self):
        percentiles = [utils.quantile(self.ranking, percentile).support() for percentile in xrange(101)]
        return percentiles

    def bgraph(self, threshold):
        bg = BGraph(self.graph, self.d, self.test_utils)
        count_true_diags = 0
        count_false_diags = 0
        count_unaligned_diags = 0
        for diag in self.__diags(threshold):
            is_true_diag = self.test_utils.is_true_diagonal(diag)
            if is_true_diag == self.test_utils.TRUE_DIAG:
              count_true_diags +=1 
            elif is_true_diag == self.test_utils.FALSE_DIAG:
              count_false_diags += 1
            else:
              count_unaligned_diags += 1
            self.test_utils.add_to_diags(diag)
            bg.add_diagonal(diag)
        self.test_utils.logger.info("true diag = "+ str(count_true_diags) + " false diag = " + str(count_false_diags) + " unanligned = " + str(count_unaligned_diags) + "\n")
        return bg
    
    def bgraph_from_genome(self):
        bg = BGraph(self.graph, self.d, self.test_utils)
        for key, rect in self.rectangles.items():
          for key, diag in rect.diagonals.items():
            bg.add_diagonal(diag)
        return bg

    def get_support(self, e1, e2):
      if (e1, e2) in self.rectangles:
        rectangle = self.rectangles[(e1,e2)]
        return max([diag.support() for diag in rectangle.diagonals.itervalues()])
      return 0

    def __diags(self, threshold = 0.0):
        assert self.ranking, "rank/filter first"
        diag_file = open("diagonals_all.txt","w")
        for d in self.ranking:
            if d.support() > threshold:
                diag_file.write(str(d.rectangle.e1.eid) + " " + str(d.rectangle.e2.eid ) + " " + str(d.D) + "\n")
                yield d
        diag_file.close()

    def __build_from_graph(self):
        for e1 in self.graph.es.itervalues():
            for e2, D in self.graph.dfs(e1, self.d):
                if (e1, e2) not in self.rectangles:
                    self.rectangles[(e1, e2)] = Rectangle(e1, e2)
                r = self.rectangles[(e1, e2)]
                r.add_diagonal(self.d, D)
                

    def __conjugate(self):
        for rect in self.rectangles.itervalues():
            conj = self.rectangles[(rect.e2.conj, rect.e1.conj)]
            conjugate(rect, conj)
            for diag in rect.diagonals.itervalues():
                assert diag.rectangle == rect
                D = diag.D - diag.rectangle.e1.len + diag.rectangle.e2.len
                pathset = diag.pathset.conj() if experimental.filter == experimental.Filter.pathsets else None
                conj = rect.conj.diagonals[D, pathset]
                conjugate(diag, conj)

    def __use_prd(self, prd_filename, config):
        for e1id, e2id, D, weight, delta in saveparser.prd(prd_filename):
            e1 = self.graph.es[e1id]
            e2 = self.graph.es[e2id]
            if (e1, e2) not in self.rectangles: # there' no possible rectangles between e1 and e2
                if weight < 0.00001 or D - e1.len + 100 < self.d:
                  continue
                if (e1,e2) not in self.not_used_prd_support:
                  self.not_used_prd_support[(e1,e2)] = []
                self.not_used_prd_support[(e1,e2)].append((D, weight,delta))
                continue
            r = self.rectangles[(e1, e2)]
            for diag in r.diagonals.itervalues():
                diag.inc_prd_support(D, weight, delta, config)
                if diag != diag.conj:
                    diag.conj.inc_prd_support(D - e1.len + e2.len, weight, delta, config) # TODO: check if prd is symmetric or not
        #print self.not_used_prd_support
     
    def use_prd_diag(self, diag):
        rect = diag.rectangle
        diag_e1_id = rect.e1.eid
        diag_e2_id = rect.e2.eid
        if (diag_e1_id, diag_e2_id) in self.prd:
           (D, weight, delta) = self.prd[(diag_e1_id, diag_e2_id)]
           diag.inc_prd_support(D, weight, delta, self.config)
           if hasattr(diag, "conj") and diag != diag.conj:
             diag.conj.inc_prd_support(D - rect.e1.len + rect.e2.len, weight, delta, self.config) 

    def __build_from_pst(self, pst_filename):
        for e1id, e2id, pathset in saveparser.pst(pst_filename):
            e1 = self.graph.es[e1id]
            e2 = self.graph.es[e2id]
            pathset = [[self.graph.es[eid] for eid in path] for path in pathset]
            self.__add_pathset(e1, e2, pathset)
        for edge in self.graph.es.itervalues():
            self.__add_pathset(edge, edge, [[edge]])

    def __add_pathset(self, e1, e2, pathset):
        if (e1, e2) not in self.rectangles:
            self.rectangles[(e1, e2)] = Rectangle(e1, e2)
        r = self.rectangles[(e1, e2)]
        # splitting pathset with different D, TODO: variance
        Ds = {} # D -> [pathset]
        for path in pathset:
            D = sum(edge.len for edge in path[:-1])
            if D in Ds:
                Ds[D] += [path]
            else:
                Ds[D] = [path]
            # making diagonal
        for D, pathset in Ds.iteritems():
            if -e2.len < D - self.d < e1.len:
                pathset = pathsets.PathSet(pathset)
                r.add_diagonal(self.d, D, pathset)

    def __rank(self):
        self.ranking = []
        for r in self.rectangles.itervalues():
            for d in r.diagonals.itervalues():
                self.ranking.append(d)
        self.ranking.sort(key=lambda x: x.support(), reverse=True)
        step = 1
        percentiles = [utils.quantile(self.ranking, q).support() for q in xrange(0,100+step,step)]
        self.logger.info("Percentiles: %s" % list(enumerate(percentiles)))

    def __all(self):
        self.ranking = [diag for r in self.rectangles.itervalues() for diag in r.diagonals.itervalues()]
