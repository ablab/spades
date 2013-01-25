############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

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
    def __init__(self, graph, d, test_utils=None, prd_file_name=None, late_pair_info_counted_prd=None, config=None):
        self.graph = graph
        self.d = d
        self.config = config
        self.prd = dict()
        self.prd_for_scaffold = dict()
        if prd_file_name:
            self.__get_prd(prd_file_name)
        if late_pair_info_counted_prd:
            self.__get_prd_for_scaffold(late_pair_info_counted_prd)
        self.rectangles = {} # (e1, e2) -> Rectangle
        self.logger = logging.getLogger('rectangles')
        self.test_utils = test_utils

    def __get_prd(self, prd_file_name):
        for e1id, e2id, D, weight, delta in saveparser.prd(prd_file_name):
            if (e1id, e2id) not in self.prd:
                self.prd[(e1id, e2id)] = []
            self.prd[(e1id, e2id)].append((D, weight, delta))

    def __get_prd_for_scaffold(self, prd_file_name):
        for e1id, e2id, D, weight, delta in saveparser.prd(prd_file_name):
            if (e1id, e2id) not in self.prd:
                e1 = self.graph.es[e1id]
                e2 = self.graph.es[e2id]
                if len(e1.v2.out) != 0 or len(e2.v1.inn) != 0:
                    continue
                if D - e1.len > 0 and D - e1.len < self.d:
                    if (e1id, e2id) not in self.prd_for_scaffold:
                        self.prd_for_scaffold[(e1id, e2id)] = []
                    self.prd_for_scaffold[(e1id, e2id)].append((D, weight, delta))

    def filter(self, prd_filename, config):
        self.__build_from_graph(config)
        #self.__conjugate()
        #self.__use_prd(prd_filename, config)
        #self.__rank()
    
    def bgraph(self, threshold):
        bg = BGraph(self.graph, self.d, self.test_utils)
        for diag in self.__diags(threshold):
            bg.add_diagonal(diag)
        
        return bg

    def __diags(self, threshold=0.0):
        #assert self.ranking, "rank/filter first"
        for d in [diag for r in self.rectangles.itervalues() for diag in r.diagonals.itervalues()]:
        #for d in self.ranking#:
            if d.support() > threshold:
                yield d

    def __build_from_graph(self, config):
        for e1 in self.graph.es.itervalues():
            for e2, D in self.graph.dfs(e1, self.d):
                if e1.eid != e2.eid and (e1.eid, e2.eid) not in self.prd:
                    continue
                if (e1, e2) not in self.rectangles:
                    self.rectangles[(e1, e2)] = Rectangle(e1, e2)
                r = self.rectangles[(e1, e2)]
                if (e2.conj, e1.conj) not in self.rectangles:
                    self.rectangles[(e2.conj, e1.conj)] = Rectangle(e2.conj, e1.conj)
                r_conj = self.rectangles[(e2.conj, e1.conj)]
                conjugate(r, r_conj)
                r.add_diagonal(self.d, D)
                D_conj = D - e1.len + e2.len
                r_conj.add_diagonal(self.d, D_conj)
                diag = r.diagonals[(D, None)]
                diag_conj = r_conj.diagonals[(D_conj, None)]
                conjugate(diag, diag_conj)
                if (e1.eid, e2.eid) in self.prd:
                    for (D1, weight, delta) in self.prd[(e1.eid, e2.eid)]:
                        diag.inc_prd_support(D1, weight, delta, config)
                        if diag != diag.conj:
                            diag.conj.inc_prd_support(D1 - e1.len + e2.len, weight, delta, config)     
                if (e2.conj.eid, e1.conj.eid) in self.prd:
                    for (D1, weight, delta) in self.prd[(e2.conj.eid, e1.conj.eid)]:
                        diag_conj.inc_prd_support(D1, weight, delta, config)
                        if diag != diag.conj:
                            diag.inc_prd_support(D1 - e2.conj.len + e1.conj.len, weight, delta, config)     
                    
                if diag.support() < 0.0000000001:
                    if (D,None) in r.diagonals:
                        del r.diagonals[(D, None)]
                    if (D_conj, None) in r_conj.diagonals:
                        del r_conj.diagonals[(D_conj, None)]
                if len(r.diagonals.keys()) == 0:
                    if (e1,e2) in self.rectangles:
                        del self.rectangles[(e1,e2)]
                    if (e2.conj, e1.conj) in self.rectangles:
                        del self.rectangles[(e2.conj, e1.conj)]



                    
    
    def __conjugate(self):
        for rect in self.rectangles.itervalues():
            conj = self.rectangles[(rect.e2.conj, rect.e1.conj)]
            conjugate(rect, conj)
            for diag in rect.diagonals.itervalues():
                assert diag.rectangle == rect
                D = diag.D - diag.rectangle.e1.len + diag.rectangle.e2.len
                if experimental.filter == experimental.Filter.pathsets:
                    pathset = diag.pathset.conj()
                else:
                    pathset = None
                conj = rect.conj.diagonals[D, pathset]
                conjugate(diag, conj)

    def __use_prd(self, prd_filename, config):
        for e1id, e2id, D, weight, delta in saveparser.prd(prd_filename):
            e1 = self.graph.es[e1id]
            e2 = self.graph.es[e2id]
            if (e1, e2) not in self.rectangles: # there' no possible rectangles between e1 and e2
                continue
            r = self.rectangles[(e1, e2)]
            for diag in r.diagonals.itervalues():
                diag.inc_prd_support(D, weight, delta, config)
                if diag != diag.conj:
                    diag.conj.inc_prd_support(D - e1.len + e2.len, weight, delta,
                        config) # TODO: check if prd is symmetric or not
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

    
    def __rank(self):
        self.ranking = []
        for r in self.rectangles.itervalues():
            for d in r.diagonals.itervalues():
                self.ranking.append(d)
        self.ranking.sort(key=lambda x: x.support(), reverse=True)
        step = 1
        percentiles = [utils.quantile(self.ranking, q).support() for q in xrange(0, 100 + step, step)]
        self.logger.info("Percentiles: %s" % list(enumerate(percentiles)))

    
    #next code about pathset, not_used
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

    def pathsets(self, pst_filename):
        self.__build_from_pst(pst_filename)
        self.__conjugate()
        self.__all()

    def __all(self):
        self.ranking = [diag for r in self.rectangles.itervalues() for diag in r.diagonals.itervalues()]
