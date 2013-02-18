############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from abstract_graph import Abstract_Graph, Abstract_Vertex, Abstract_Edge
import logging
import N50
import sys
import utils
import saveparser
from utils import conjugate

##################
# DEBRUIJN GRAPH #
##################

class Vertex(Abstract_Vertex):
    def __init__(self, vid, conj):
        Abstract_Vertex.__init__(self, vid)
        conjugate(self, conj)

    def seq(self, K):
        if self.out:
            return self.out[0].seq[:K]
        else:
            return self.inn[0].seq[-K:]

    def __repr__(self):
        return "V%d" % (self.vid)


class Edge(Abstract_Edge):
    def __init__(self, eid, v1, v2, edge_len, conj):
        Abstract_Edge.__init__(self, v1, v2, eid)
        self.len = edge_len
        conjugate(self, conj)
        self.seq = None
        self.cvr = 0

    def length(self):
        return self.len

    def __repr__(self):
        return "E%d(%d)" % (self.eid, self.len)
    
    def print_fasta(self, stream, contig_id_start):
        contig_id = contig_id_start
        for id_contig, contig in enumerate(self.seq.split('N')):
            if not contig: 
                continue
            l = len(contig)
            print >> stream, '>contig_%d_%d_%d_%d_l=%06d' % (contig_id, self.eid, self.conj.eid, id_contig, l)
            contig_id += 1
            for l in xrange(0, l, 60):
                print >> stream, contig[l:l + 60]
            
        return contig_id

       
        

class Graph(Abstract_Graph):
    def __init__(self):
        Abstract_Graph.__init__(self)
        self.max_eid = 0
        self.etalon_dist = dict()
        self.logger = logging.getLogger('rectangles')

    def add_vertex(self, vid, conj_id):
        #assert vid != conj_id, "Vertex can't be self-conjugated"
        conj = self.vs.get(conj_id, None)
        v = Vertex(vid, conj)
        self.vs[vid] = v
        return v

    def add_edge(self, eid, v1id, v2id, edge_len, conj_id):
        #assert eid != conj_id, "Self-conjugate edges are not supported yet"
        if eid in self.es:
            return self.es[eid]
        if eid > self.max_eid or conj_id > self.max_eid:
            self.max_eid = max(eid, conj_id)
        conj = self.es.get(conj_id, None)
        v1 = self.vs[v1id]
        v2 = self.vs[v2id]
        e = Edge(eid, v1, v2, edge_len, conj)
        self.es[eid] = e
        if eid == conj_id:
            conjugate(e, e)
        return e

    def add_seq(self, eid, seq):
        self.es[eid].seq = seq

    def add_cvr(self, eid, cvr):
        self.es[eid].cvr = cvr

    def update_K(self):
        assert len(self.es) > 0, "Empty graph"
        any_edge = (e for e in self.es.itervalues()).next()
        K = len(any_edge.seq) - any_edge.len
        self.K = K

    def check(self):
        for v in self.vs.itervalues():
            assert v.conj, "Some vertex have no conjugate"
        for e in self.es.itervalues():
            assert e.conj, "Some edge have no conjugate"
        for e in self.es.itervalues():
            assert self.K == len(e.seq) - e.len, "Inconsistent K"
        for e in self.es.itervalues():
            assert e.seq == utils.rc(e.conj.seq), (e.seq, utils.rc(e.conj.seq))

    def find_loops(self, threshold, L, rs):
        edges_before_loop = dict()
        for eid, e in self.es.items():
            if e.len > L:
                result_loop = self.find_all_loops(e, threshold, L, rs.rectangles)
                if result_loop:
                    edges_before_loop[e.eid] = result_loop
        return edges_before_loop

    def fasta_for_long_contigs(self, K, d, is_sc, is_careful,  stream=sys.stdout, should_connect=dict(), scaffold=dict()):
        in_paths = set()
        for edge_id, path in should_connect.items():
            in_paths.add(path[-1].eid)
            in_paths.add(path[-1].conj.eid)
            in_paths.add(path[0].eid)
            in_paths.add(path[0].conj.eid)
            """for e in path:
            in_paths.add(e.eid)
            in_paths.add(e.conj.eid)"""
        for edge_id, path in scaffold.items():
            for e in path:
                in_paths.add(e.eid)
                in_paths.add(e.conj.eid)
        contig_id = 0
        
        empty_format_function = lambda x : ""
        scaffold_format_function = lambda be : be.get_seq_for_contig(K, d, is_sc, is_careful) + "NNN"
        should_connect_pre_format = lambda path : path[0].get_begin_seq(K, d, is_sc)
        should_connect_format = lambda be: be.get_midle_seq() 
        should_connect_post_format = lambda path : path[-1].get_end_seq(K, d, is_sc)

        for edge in self.es.itervalues():
            if edge.conj.eid <= edge.eid: # non-conjugate
                if (print_graph_edge(edge, scaffold, stream, empty_format_function, scaffold_format_function, empty_format_function) or
                        print_graph_edge(edge, should_connect, stream, should_connect_pre_format, should_connect_format, should_connect_post_format)):
                    continue
                
                if edge.eid in in_paths:
                    continue
                
                contig_id = edge.print_fasta(stream, contig_id)

    def fasta(self, stream=sys.stdout):
        contig_id = 0
        for edge in self.es.itervalues():
            if edge.conj.eid <= edge.eid: # non-conjugate
                contig_id = edge.print_fasta(stream, contig_id)
    
    
    def stats(self, d):
        ls = []
        for edge in self.es.itervalues():
            if edge.conj.eid <= edge.eid: # non-conjugate
                ls.append(len(edge.seq))
        ls.sort()
        self.logger.info('Edges   = %s' % ls)
        self.logger.info('#Edges  = %d' % len(ls))
        self.logger.info('Edg N50 = %dbp' % N50.N50(ls))
        ls = []
        for edge in self.es.itervalues():
            #if 'N' in edge.seq: continue
            if edge.conj.eid <= edge.eid: # non-conjugate
                lens = filter(None, map(len, edge.seq.split('N')))
                ls += lens
        ls.sort()
        self.logger.info('Contigs = %s' % ls)
        self.logger.info('#Contigs= %d' % len(ls))
        self.logger.info('K       = %d' % self.K)
        self.logger.info('d       = %d' % d)
        self.logger.info('Total   = %dbp' % sum(ls))
        self.logger.info('N50     = %d' % N50.N50(ls))
        self.logger.info('We split small edges (with Ns) to multiple contigs, so #edges < #contigs')
        self.logger.info('N50>1000 = %d' % N50.N50([x for x in ls if x >= 1000]))
        return N50.N50(ls)

    def save(self, filename):
        # Graph save
        grp = open(filename + '.grp', 'w')
        print >> grp, len(self.vs), len(self.es)
        for vertex in self.vs.itervalues(): # Vertex 2 ~ 1 .
            print >> grp, 'Vertex %d ~ %d .' % (vertex.vid, vertex.conj.vid)#, vertex.inn, vertex.out, vertex.seq(15)
        print >> grp  # empty line
        for edge in self.es.itervalues():  # Edge 15 : 12 -> 2, l = 838 ~ 16 .
            if not edge or not edge.conj:
                continue
            print >> grp, 'Edge %d : %d -> %d, l = %d ~ %d .' % (
                edge.eid, edge.v1.vid, edge.v2.vid, len(edge.seq), edge.conj.eid)#, self.etalon_dist[edge.eid], edge.seq
        grp.close()

        # Sequences save
        sqn = open(filename + '.sqn', 'w')
        contigs = open(filename + ".contigs", "w")
        print >> sqn, len(self.es)
        visited = set()
        for edge in self.es.itervalues():  # Edge 15 : 12 -> 2, l = 838 ~ 16 .
            if edge.eid not in visited:
                visited.add(edge.eid)
                visited.add(edge.conj.eid)
                print >> contigs, ">" + str(edge.eid) + "\n" + edge.seq.strip()
            print >> sqn, '%d %s .' % (edge.eid, edge.seq)
        sqn.close()
        contigs.close()
        # Coverage save
        cvr = open(filename + '.cvr', 'w')
        print >> cvr, len(self.es)
        for edge in self.es.itervalues():  # 15 1.234 .
            print >> cvr, '%d %f .' % (edge.eid, edge.cvr)
        cvr.close()
        # Stupid pos file
        pos = open(filename + '.pos', 'w')
        print >> pos, len(self.es)
        for edge in self.es.itervalues():  # 15 0
            print >> pos, '%d %d' % (edge.eid, 0)
        pos.close()

    def load(self, grp_filename, sqn_filename, cvr_filename):
        for vid, conj in saveparser.grp_vertices(grp_filename):
            self.add_vertex(vid, conj)
        for eid, v1id, v2id, l, conj in saveparser.grp_edges(grp_filename):
            self.add_edge(eid, v1id, v2id, l, conj)
        for eid, seq in saveparser.sqn(sqn_filename):
            self.add_seq(eid, seq)
        for eid, cvr in saveparser.cvr(cvr_filename):
            self.add_cvr(eid, cvr)
        self.update_K()

    def __get_kmers_pos(self, genome, k):
        kmers = dict()
        for i in range(len(genome) - k + 1):
            kmer = genome[i: i + k]
            if kmer not in kmers:
                kmers[kmer] = [i + 1]
                kmers[utils.rc(kmer)] = [-(len(genome) - i - k + 1)]
            else:
                kmers[kmer].append(i + 1)
                kmers[utils.rc(kmer)].append(-(len(genome) - i - k + 1))
        return kmers

    def dfs(self, e, d):
        limit1 = d - e.len
        limit2 = d
        if e.len > d:
            yield e, 0
        ls = [set() for _ in xrange(limit2)]
        ls[0].add(e.v2)
        all_dist = dict()
        for pos in xrange(limit2):
            for v in ls[pos]:
                for e2 in v.out:
                    pos2 = pos + e2.len
                    if pos2 < limit2:
                        ls[pos2].add(e2.v2)
                    if pos + e2.len > limit1:
                        yield e2, pos + e.len

def print_graph_edge(edge, paths, stream, pre_format_function, format_function, post_format_function):
    if not edge.eid in paths:
        return False    
    str_id = ""
    #print "PRINTING SCAFFOLD", edge.eid, [e.eid for e in scaffold[edge.eid]]
    path = paths[edge.eid]
    seq = pre_format_function(path)
    for be in path:
        seq += format_function(be)
        str_id += "_" + str(be.eid) + "_"
    seq += post_format_function(path)
    print >> stream, '>contig_%d_l=%06d_long%s' % (edge.eid, len(seq), str_id)
    l = len(seq)
    for l in xrange(0, l, 60):
        print >> stream, seq[l:l + 60]
    return True


