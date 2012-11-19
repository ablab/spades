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
      return self.out[0].seq[ :K]
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


class Graph(Abstract_Graph):
  def __init__(self):
    Abstract_Graph.__init__(self)
    self.max_eid = 0
    self.etalon_dist = dict()
    self.logger = logging.getLogger('rectangles')

  def add_vertex(self, vid, conj_id):
    assert vid != conj_id, "Vertex can't be self-conjugated"
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

  def find_loops(self, threshold, L):
    edges_before_loop = dict()
    for eid, e in self.es.items():
      if e.len > L:
        result_loop = self.find_all_loops(e, threshold, L)
        if result_loop:
          edges_before_loop[e.eid] = result_loop
    return edges_before_loop
  
  def fasta_for_long_contigs(self, K, d, is_sc, stream=sys.stdout, should_connect = dict(), scaffold = dict()):
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
        for edge in self.es.itervalues():
            if edge.conj.eid <= edge.eid: # non-conjugate
              if edge.eid in scaffold:
                str_id = ""
                print "PRINTING SCAFFOLD", edge.eid, [e.eid for e in scaffold[edge.eid]]
                path = scaffold[edge.eid]
                seq = ""
                for be in path:
                  seq += be.get_seq_for_contig(K, d, is_sc) + "NNN"
                  str_id += "_" + str(be.eid) + "_" 
                print >>stream,  '>contig_%d_l=%06d_long%s' % (edge.eid, len(seq), str_id)
                l = len(seq)
                for l in xrange(0, l, 60):
                    print >>stream, seq[l:l + 60]
                continue
              if edge.eid in should_connect:
                str_id = ""
                print "PRINTING", edge.eid, [e.eid for e in should_connect[edge.eid]]
                path = should_connect[edge.eid]
                seq = path[0].get_begin_seq(K, d, is_sc)
                seq = ""
                for be in path:
                  seq += be.get_midle_seq()
                  str_id += "_" + str(be.eid) + "_" 
                seq += path[-1].get_end_seq(K, d, is_sc)
                print >>stream,  '>contig_%d_l=%06d_long%s' % (edge.eid, len(seq), str_id)
                l = len(seq)
                for l in xrange(0, l, 60):
                    print >>stream, seq[l:l + 60]
                continue

              if edge.eid in in_paths:
                continue     
              for id_contig, contig in enumerate(edge.seq.split('N')):
                    if not contig: continue
                    l = len(contig)
                    print >>stream, '>contig_%d_%d_%d_%d_l=%06d' % (contig_id, edge.eid, edge.conj.eid, id_contig, l)
                    contig_id += 1
                    for l in xrange(0, l, 60):
                        print >>stream, contig[l:l + 60]

  def fasta(self, stream=sys.stdout):
    contig_id = 0 
    for edge in self.es.itervalues():
      if edge.conj.eid <= edge.eid: # non-conjugate
        for id_contig, contig in enumerate(edge.seq.split('N')):
          if not contig: continue
          l = len(contig)
          print >>stream, '>contig_%d_%d_%d_%d_l=%06d' % (contig_id, edge.eid, edge.conj.eid, id_contig, l)
          contig_id += 1
          for l in xrange(0, l, 60):
            print >>stream, contig[l:l + 60]

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
    print >>grp, len(self.vs), len(self.es)
    for vertex in self.vs.itervalues(): # Vertex 2 ~ 1 .
      print >>grp, 'Vertex %d ~ %d .' % (vertex.vid, vertex.conj.vid)#, vertex.inn, vertex.out, vertex.seq(15)
    print >>grp  # empty line
    for edge in self.es.itervalues():  # Edge 15 : 12 -> 2, l = 838 ~ 16 .
      if not edge or not edge.conj:
        continue
      print >>grp, 'Edge %d : %d -> %d, l = %d ~ %d .' % (
        edge.eid, edge.v1.vid, edge.v2.vid, len(edge.seq), edge.conj.eid)#, self.etalon_dist[edge.eid], edge.seq
    grp.close()

    # Sequences save
    sqn = open(filename + '.sqn', 'w')
    contigs = open(filename + ".contigs", "w")
    print >>sqn, len(self.es)
    visited = set()
    for edge in self.es.itervalues():  # Edge 15 : 12 -> 2, l = 838 ~ 16 .
      if edge.eid not in visited:
        visited.add(edge.eid)
        visited.add(edge.conj.eid)
        print >> contigs, ">" + str(edge.eid) + "\n" + edge.seq.strip()
      print >>sqn, '%d %s .' % (edge.eid, edge.seq)
    sqn.close()
    contigs.close()
    # Coverage save
    cvr = open(filename + '.cvr', 'w')
    print >>cvr, len(self.es)
    for edge in self.es.itervalues():  # 15 1.234 .
      print >>cvr, '%d %f .' % (edge.eid, edge.cvr)
    cvr.close()
    # Stupid pos file
    pos = open(filename + '.pos', 'w')
    print >>pos, len(self.es)
    for edge in self.es.itervalues():  # 15 0
      print >>pos, '%d %d' % (edge.eid, 0)
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
      kmer = genome[i: i +k]
      if kmer not in kmers:
        kmers[kmer] = [i+1]
        kmers[utils.rc(kmer)] = [-(len(genome) - i -k + 1)]
      else:
        kmers[kmer].append(i+1)
        kmers[utils.rc(kmer)].append(-(len(genome)- i - k+1))
    return kmers

  def make_graph(self, genome, k):
    self.K = k
    kmers = self.__get_kmers_pos(genome, k)
    visit = set()
    vid = 0
    eid = 0
    edges = set()
    verts = dict()
    for key in kmers:
      if key in visit:
        continue
      body = [key[-1]]
      end_vertex = key[1:]                    
      while True:
          next_kmer = extend_forward(end_vertex, kmers)
          if next_kmer == None:
              break
          body.append(next_kmer[-1])
          end_vertex = next_kmer[1:]
          visit.add(next_kmer)
          visit.add(utils.rc(next_kmer))
          
      begin_vertex = key[:-1]
      while True:
          next_kmer = extend_backward(begin_vertex, kmers)
          if next_kmer == None:
              break
          body.insert(0, next_kmer[-1])
          begin_vertex = next_kmer[0:-1]
          visit.add(next_kmer)
          visit.add(utils.rc(next_kmer))
          
      body = begin_vertex + ''.join(body)
      if begin_vertex not in verts:
        begin_ref = self.add_vertex(vid, vid+1)
        r_end_ref = self.add_vertex(vid+1, vid)
        verts[begin_vertex] = begin_ref.vid
        verts[utils.rc(begin_vertex)] = r_end_ref.vid
        vid +=2
      if end_vertex not in verts:
        end_ref = self.add_vertex(vid, vid+1)
        r_begin_ref = self.add_vertex(vid+1, vid)
        verts[end_vertex] = end_ref.vid
        verts[utils.rc(end_vertex)] = r_begin_ref.vid
        vid +=2
      bv = verts[begin_vertex]
      ev = verts[end_vertex]
      rbv = verts[utils.rc(end_vertex)]
      rev = verts[utils.rc(begin_vertex)]
      if (bv, ev) not in edges:
        if (bv,ev) == (rbv, rev) and body == utils.rc(body):
          self.add_edge(eid, bv, ev, len(body) -k +1 , eid)
          edges.add((bv,ev))
          self.add_seq(eid, body)
          self.etalon_dist[eid] = kmers[body[:k]] + kmers[utils.rc(body)[:k]]
          eid += 1
        else:
          self.add_edge(eid, bv, ev, len(body) - k + 1, eid +1)
          self.add_edge(eid +1, rbv, rev, len(body) -k +1, eid)
          edges.add((bv,ev)) 
          edges.add((rbv, rev))
          self.add_seq(eid, body)
          self.add_seq(eid +1, utils.rc(body))
          self.etalon_dist[eid] = kmers[body[:k]]
          self.etalon_dist[eid+1] = kmers[utils.rc(body)[:k]]
          eid += 2
  
  def __from_genome(self):
    return len(self.etalon_dist.keys()) > 0

  def dfs(self, e, d):
    limit1 = d - e.len
    limit2 = d
    if e.len > d:
      yield e, 0
    ls = [set() for _ in xrange(limit2)]
    ls[0].add(e.v2)
    all_dist  = dict()
    if self.__from_genome():
      all_dist[(0,e.v2.vid)] = (e, self.etalon_dist[e.eid])
    for pos in xrange(limit2):
      for v in ls[pos]:
        if self.__from_genome():
          (prev_e, dists) = all_dist[(pos, v.vid)]
        for e2 in v.out:
          if self.__from_genome():
            new_dists = []
            for dist in dists:
              if dist>= 0 and dist + prev_e.len  + 1 in self.etalon_dist[e2.eid]:
                new_dists.append(dist + prev_e.len + 1)
              elif dist<= 0 and -(-dist + prev_e.len  + 1) in self.etalon_dist[e2.eid]:
                new_dists.append(dist - prev_e.len - 1)
            if len(new_dists) == 0:
              continue
          pos2 = pos + e2.len
          if pos2 < limit2:
            ls[pos2].add(e2.v2)
            if self.__from_genome():
              if (pos2, e2.v2.vid) in all_dist:
                new_dists += all_dist[(pos2, e2.v2.vid)][1]
              all_dist[(pos2, e2.v2.vid)] = (e2, new_dists)
          if pos + e2.len > limit1:
            yield e2, pos + e.len

alphabet = "ACGT"

def find_out_edges(vertex_body, kmer_map):
  next_kmer = [(vertex_body + base) for base in alphabet]
  return [kmer for kmer in next_kmer if kmer in kmer_map]

def find_in_edges(vertex_body, kmer_map):
  next_kmer = [(base + vertex_body) for base in alphabet]
  return [kmer for kmer in next_kmer if kmer in kmer_map]

def extend_forward(vertex_body, kmer_map):    
  in_edge = find_in_edges(vertex_body, kmer_map)
  out_edge = find_out_edges(vertex_body, kmer_map)
  if len(in_edge) == 1 and len(out_edge) == 1:
    return out_edge[0]
  return None

def extend_backward(vertex_body, kmer_map):    
  in_edge = find_in_edges(vertex_body, kmer_map)
  out_edge = find_out_edges(vertex_body, kmer_map)
  if len(in_edge) == 1 and len(out_edge) == 1:
    return in_edge[0]
  return None

    
