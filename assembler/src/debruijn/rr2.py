#!/usr/bin/python

###########################################################################################
###
### SPAdes Repeat Resolver based on PNAS paper
###
### INPUT: De Bruijn Graph + Paired Info + Edge Sequences
### OUTPUT: FASTA, some Graph
###
### TODO: what type of contigs' prolongation should we use? do we need scaffolding?
### TODO: check outputted graph
### TODO: check misassemblies on SC_LANE1
###
###########################################################################################

import sys
sys.path.append('src/tools/quality/libs') # for N50
import N50
import parser
import itertools
import cStringIO

##################
# UTIL FUNCTIONS #
##################

def conjugate(a, b):
    a.conj = b
    if b: b.conj = a

def seq_join(seq1, seq2): # join two sequences, prioritize ACGT over N
    assert len(seq1) == len(seq2)
    for c1, c2 in zip(seq1, seq2):
        assert c1 == c2 or c1 == 'N' or c2 == 'N', seq1 + '|' + seq2
    seq = ''.join(map(lambda x, y: x if x != 'N' else (y if y != 'N' else 'N'), seq1, seq2))
    return seq

##################
# DEBRUIJN GRAPH #
##################

class Graph(object):
    class Vertex(object):
        def __init__(self, vid, conj):
            self.vid = vid
            self.inn = []
            self.out = []
            conjugate(self, conj)

        def __hash__(self):
            return self.vid

        def seq(self, K):
            return self.out[0].seq[ :K]

    class Edge(object):
        def __init__(self, eid, v1, v2, len, conj):
            self.eid = eid
            self.v1 = v1
            self.v2 = v2
            self.len = len
            conjugate(self, conj)
            self.seq = None

        def __hash__(self):
            return self.eid

    def __init__(self):
        self.vs = {} # vid -> Vertex
        self.es = {} # eid -> Edge
        self.max_eid = 0

    def add_vertex(self, vid, conj_id):
        assert vid != conj_id, "Vertex can't be self-conjugated"
        conj = self.vs.get(conj_id, None)
        v = Graph.Vertex(vid, conj)
        self.vs[vid] = v
        return v

    def add_edge(self, eid, v1id, v2id, len, conj_id):
        #assert eid != conj_id, "Self-conjugate edges are not supported yet"
        if eid > self.max_eid or conj_id > self.max_eid:
            self.max_eid = max(eid, conj_id)
        conj = self.es.get(conj_id, None)
        v1 = self.vs[v1id]
        v2 = self.vs[v2id]
        e = Graph.Edge(eid, v1, v2, len, conj)
        v1.out.append(e)
        v2.inn.append(e)
        self.es[eid] = e
        if eid == conj_id:
            conjugate(e, e)
        return e

    def add_N_edge(self, v1id, v2id, len, K, conj_id):
        self.max_eid += 1
        e = self.add_edge(self.max_eid, v1id, v2id, len, conj_id)
        seq1 = self.vs[v1id].seq(K) + 'N'*len
        seq2 = 'N'*len + self.vs[v2id].seq(K)
        seq = seq_join(seq1, seq2)
        self.add_seq(self.max_eid, seq)
        return e

    def add_two_N_edges(self, v1, v2, len, K):
        e1 = self.add_N_edge(v1.vid, v2.vid, len, K, None)
        e2 = self.add_N_edge(v2.conj.vid, v1.conj.vid, len, K, e1.eid)
        return e1, e2

    def add_seq(self, eid, seq):
        self.es[eid].seq = seq

    def get_edge(self, eid):
        return self.es[eid]

    def update_K(self):
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

    def fasta(self, stream=sys.stdout):
        for edge in self.es.itervalues():
            if edge.conj.eid <= edge.eid: # non-conjugate
                l = len(edge.seq)
                print >> stream, '>contig_%d_%d_l=%06d' % (edge.eid, edge.conj.eid, l)
                for l in xrange(0, l, 60):
                    print >> stream, edge.seq[l:l + 60]

    def stats(self, stream=sys.stderr):
        ls = []
        for edge in self.es.itervalues():
            if edge.conj.eid <= edge.eid: # non-conjugate
                ls.append(len(edge.seq))
        ls.sort()
        print >> stream, 'Ls  =', ls
        print >> stream, 'K   =', self.K
        print >> stream, 'Len =', sum(ls), 'bp'
        print >> stream, 'N50 =', N50.N50(ls), 'bp'
        print >> stream, 'Num =', len(ls), 'contigs'

    def save(self, filename):
        # Graph save
        grp = open(filename + '.grp', 'w')
        print >> grp, len(self.vs), len(self.es)
        for vertex in self.vs.itervalues(): # Vertex 2 ~ 1 .
            print >> grp, 'Vertex %d ~ %d .' % (vertex.vid, vertex.conj.vid)
        print >> grp  # empty line
        for edge in self.es.itervalues():  # Edge 15 : 12 -> 2, l = 838 ~ 16 .
            print >> grp, 'Edge %d : %d -> %d, l = %d ~ %d .' % (
                edge.eid, edge.v1.vid, edge.v2.vid, len(edge.seq), edge.conj.eid)
        grp.close()
        # Sequences save
        sqn = open(filename + '.sqn', 'w')
        print >> sqn, len(self.es)
        for edge in self.es.itervalues():  # Edge 15 : 12 -> 2, l = 838 ~ 16 .
            print >> sqn, '%d %s .' % (edge.eid, edge.seq)
        sqn.close()

    def load(self, grp_filename, sqn_filename):
        for vid, conj in parser.grp_vertices(grp_filename):
            self.add_vertex(vid, conj)
        for eid, v1id, v2id, l, conj in parser.grp_edges(grp_filename):
            self.add_edge(eid, v1id, v2id, l, conj)
        for eid, seq in parser.sqn(sqn_filename):
            self.add_seq(eid, seq)
        self.update_K()

#################
#    DIAGONAL   #
#################

class Diagonal(object):

    def __init__(self, e1, e2, D):
        # Diagonal from (e1, offseta | e2, offsetb) to (e1, offsetc | e2, offsetd)
        self.e1, self.e2, self.D = e1, e2, D
        self.offseta, self.offsetb, self.offsetc, self.offsetd = Diagonal.offsets(e1, e2, D)
        self.key1 = Diagonal.key(self.e1, self.offseta, self.e2, self.offsetb)
        self.key2 = Diagonal.key(self.e1, self.offsetc, self.e2, self.offsetd)

    @staticmethod
    def offsets(e1, e2, D):
        l1, l2 = e1.len, e2.len
        # calculate offsets in rectangle
        if d >= D:
            offseta, offsetb = 0, d - D
        else:
            offseta, offsetb = D - d, 0
        if d >= D + l2 - l1:
            offsetc, offsetd = l2 + (D - d), l2
        else:
            offsetc, offsetd = l1, l1 + (d - D)
        assert offsetc - offseta == offsetd - offsetb, "Should be main diagonal (it's a bug)"
        return offseta, offsetb, offsetc, offsetd

    @staticmethod
    def valid(e1, e2, D):
        # valid rectangle && not 0-length rectangle
        l1, l2 = e1.len, e2.len
        offseta, offsetb, offsetc, offsetd = Diagonal.offsets(e1, e2, D)
        return (0 <= offseta <= l1) and (0 <= offsetb <= l2) and (0 <= offsetc <= l1) and (0 <= offsetd <= l2) and (
            offsetc != offseta) and (offsetb != offsetd)

    @staticmethod
    def valid_or_angle(e1, e2, D):
        # valid rectangle && not 0-length rectangle
        l1, l2 = e1.len, e2.len
        offseta, offsetb, offsetc, offsetd = Diagonal.offsets(e1, e2, D)
        return (0 <= offseta <= l1) and (0 <= offsetb <= l2) and (0 <= offsetc <= l1) and (0 <= offsetd <= l2)

    @staticmethod
    def key(e1, offset1, e2, offset2):
        # key = ((e, o), v) or (v, (e, o)) or (v, v)
        first = e1.v1 if offset1 == 0 else ( e1.v2 if offset1 == e1.len else (e1, offset1) )
        second = e2.v1 if offset2 == 0 else ( e2.v2 if offset2 == e2.len else (e2, offset2) )
        return first, second

##################
#    RECTANGLE   #
##################

class Rectangle(object):
    rid = 0

    def __init__(self, e1, e2, Ds):
        # Rectangle between e1 and e2 with distances Ds
        self.e1, self.e2, self.Ds = e1, e2, Ds
        self.diagonals = []
        for D in self.Ds:
            if Diagonal.valid(e1, e2, D):
                self.diagonals.append(Diagonal(e1, e2, D))
        self.conj = None
        self.update_keys()
        self.rid = Rectangle.rid
        Rectangle.rid += 1

    def __hash__(self):
        return self.rid

    def get_keys1(self): # down-left vertex
        keys = set()
        for diagonal in self.diagonals:
            keys.add(diagonal.key1)
        return keys

    def get_keys2(self): # up-right vertex
        keys = set()
        for diagonal in self.diagonals:
            keys.add(diagonal.key2)
        return keys

    def update_keys(self):
        self.keys1, self.keys2 = self.get_keys1(), self.get_keys2()

    @staticmethod
    def valid(e1, e2, Ds):
#        for D in Ds:
#            if Diagonal.valid(e1, e2, D):
#                return True
#        return False
        for D in Ds:
            if not Diagonal.valid(e1, e2, D):
                return False
        return True

##################
#    BI GRAPH    #
##################

class BGraph(object):

    class BVertex(object):
        bvid = 0

        def __init__(self, keys):
            self.initial_keys = keys
            self.inn, self.out = [], []
            self.conj = None
            self.bvid = BGraph.BVertex.bvid
            BGraph.BVertex.bvid += 1

        def __hash__(self):
            return self.bvid

        def make_set(self): # for union-find, for gluing
            self.parent = self

        def find(self): # for union-find, for gluing
            if self.parent != self:
                self.parent = self.parent.find()
            return self.parent

        def union(self, other): # glue, for union-find
            other.find().parent = self.find()
            #self.find().initial_keys |= other.find().initial_keys # TODO: check correctness, seems we don't need it

    class BEdge(object):
        def __init__(self, rectangle, bv1, bv2):
            self.rectangles = [rectangle]
            self.beid = rectangle.rid
            self.bv1, self.bv2 = bv1, bv2
            self.conj = None

        def __hash__(self):
            return self.beid

        def get_seq(self, K):
            seq1 = cStringIO.StringIO()
            seq2 = cStringIO.StringIO()
            seq2.write('N' * d)
            for this, next in zip(self.rectangles, self.rectangles[1:]):
                assert this.diagonals[0].key2 == next.diagonals[0].key1, "Diagonals adjustment isn't implemented yet"
                #assert len(this.diagonals) == len(next.diagonals)
                diagonal = this.diagonals[0]
                # assert this.e1 != next.e1 or this.e2 != next.e2 # can happens on single-cell :)
                assert len(this.e1.seq) >= diagonal.offsetc
                assert len(this.e2.seq) >= diagonal.offsetd
                assert diagonal.offsetc - diagonal.offseta == diagonal.offsetd - diagonal.offsetb
                seq1.write(this.e1.seq[diagonal.offseta : diagonal.offsetc])
                seq2.write(this.e2.seq[diagonal.offsetb : diagonal.offsetd])
            last = self.rectangles[-1]
            diagonal = last.diagonals[0]
            seq1.write(last.e1.seq[diagonal.offseta : diagonal.offsetc + K])
            seq2.write(last.e2.seq[diagonal.offsetb : diagonal.offsetd + K])
            seq1.write('N' * d)
            #seq = seq1.getvalue().strip('N')
            seq = seq_join(seq1.getvalue(), seq2.getvalue()).strip('N')
            #assert len(seq1.getvalue()) == len(seq2.getvalue()) >= len(seq)
            return seq

    def __init__(self, graph):
        self.graph = graph
        self.rectangles = {} # (e1, e2, Ds) -> Rectangle
        self.vs = set() # BVertices
        self.es = set() # BEdges

    def add_rectangle(self, e1, e2, Ds):
        rect_Dkey = tuple(sorted(Ds))
        conj_Dkey = tuple(sorted(x - e1.len + e2.len for x in Ds))
        rect_key = (e1, e2, rect_Dkey)
        conj_key = (e2.conj, e1.conj, conj_Dkey)
        if (rect_key in self.rectangles) or not Rectangle.valid(*rect_key):
            return
        assert (conj_key not in self.rectangles) and Rectangle.valid(*conj_key), "Bug"
        rect = Rectangle(*rect_key)
        self.rectangles[rect_key] = rect
        if rect_key == conj_key:
            conj = rect
        else:
            conj = Rectangle(*conj_key)
            self.rectangles[conj_key] = conj
        conjugate(rect, conj)
        return rect, conj

    def load(self, prd_filename): # load paired info from file to rectangles
        for e1id, e2id, Ds, weight in parser.prd(prd_filename):
            e1 = self.graph.get_edge(e1id)
            e2 = self.graph.get_edge(e2id)
            #assert (e1 != e2.conj) and (e2 != e1.conj), "Self-conjugate rectangles aren't supported (even partly)"
            self.add_rectangle(e1, e2, Ds)
            # manually add all self-rectangles (e-to-e), just to be sure (if they are missing in prd file)
        for e in self.graph.es.itervalues():
            self.add_rectangle(e, e, [0])

    def build(self): # build simple bi-edges and bi-vertices from rectangles
        map = {} # rectangle -> BiEdge, for conjugativity
        for rectangle in self.rectangles.itervalues():
            bv1 = BGraph.BVertex(rectangle.keys1)
            bv2 = BGraph.BVertex(rectangle.keys2)
            be = BGraph.BEdge(rectangle, bv1, bv2)
            self.vs.add(bv1)
            self.vs.add(bv2)
            map[rectangle] = be
        for bedge in map.itervalues():
            assert len(bedge.rectangles) == 1
            conj = map[bedge.rectangles[0].conj]
            conjugate(bedge, conj)
            conjugate(bedge.bv1, conj.bv2)
            conjugate(bedge.bv2, conj.bv1)
        self.es = set(map.itervalues())

    def scaffold(self): # scaffold two tips by an edge
        scaffolds = 0
        new_bes = set()
        for be1 in bgraph.es:
            assert len(be1.rectangles) == 1, "Should be run before condensation"
            if be1.bv2.out: continue
            for be2 in bgraph.es:
                assert len(be2.rectangles) == 1, "Should be run before condensation"
                if be2.bv1.inn: continue
                if be2 == be1: continue
                # ---be1(r1)---> bv1 , bv2 ---be2(r2)--->
                r1 = be1.rectangles[0]
                r2 = be2.rectangles[0]
                bv1 = be1.bv2
                bv2 = be2.bv1
                if r1.e1 != r2.e1: continue
                # transforms to:
                # ---be1(r1)---> bv1 ----new_be1(new_r1)--->  bv2 ---be2(r2)--->
                #                    <---new_be2(new_r2)----
                if len(r1.diagonals) != len(r2.diagonals): continue # don't know how to scaffold two rectangles with different number of diagonals
                diffs = map(lambda diag1, diag2: diag2.offseta - diag1.offsetc, r1.diagonals, r2.diagonals)
                diff = diffs[0]
                if diffs.count(diff) != len(diffs): continue # don't know how to scaffold two rectangles with different number of diffs
                if diff < 0 or diff > max_scaffold: continue
                assert diff != 0, "Strange, they should be already glued"
                new_e1, new_e2 = self.graph.add_two_N_edges(r1.e2.v2, r2.e2.v1, diff, self.graph.K)
                new_r1, new_r2 = self.add_rectangle(r1.e1, new_e1, [d + diag.offsetc for diag in r1.diagonals])
                new_be1 = BGraph.BEdge(new_r1, bv1, bv2)
                new_be2 = BGraph.BEdge(new_r2, bv2.conj, bv1.conj)
                conjugate(new_be1, new_be2)
                new_bes.add(new_be1)
                new_bes.add(new_be2)
                bv1.out.append(new_be1)
                bv2.inn.append(new_be1)
                bv2.conj.out.append(new_be2)
                bv1.conj.inn.append(new_be2)
                scaffolds += 1
        bgraph.es |= new_bes
        print >>logfile, "Scaffolding: %d pairs of rectanges were added during scaffolding" % scaffolds

    def glue(self): # glue bi-vertices based on keys intersection
        map = {} # key -> BVertex
        for bvertex in self.vs:
            bvertex.make_set()
        for bvertex in self.vs:
            for key in bvertex.initial_keys:
                if key in map and map[key].find() != bvertex.find():
                    bvertex.union(map[key]) # glue them (union)!
                    bvertex.conj.union(map[key].conj)
                map[key] = bvertex.find()
        for bedge in self.es:
            bedge.bv1 = bedge.bv1.find()
            bedge.bv2 = bedge.bv2.find()
            bedge.bv1.out.append(bedge)
            bedge.bv2.inn.append(bedge)
        l = len(self.vs)
        self.vs = set(bv for bv in self.vs if len(bv.inn) > 0 or len(bv.out) > 0)
        print >>logfile, "Glued %d bi-vertices (out of %d). %d bi-vertices left." % (
            l - len(self.vs), l, len(self.vs))

    def join_biedges(self, be1, be2):
        assert be1.bv2 == be2.bv1
        be3 = be2.conj
        be4 = be1.conj
        u, v, w = be1.bv1, be1.bv2, be2.bv2
        x, y, z = be3.bv1, be3.bv2, be4.bv2
        assert 1 == len(v.inn) == len(v.out) == len(y.inn) == len(y.out)
        ## u ---be1---> v ---be2---> w
        ## z <--be4---- y <--be3---- x
        ## transforms to:
        ## u --------be1--------> w
        ## z <-------be3--------- x
        be1.rectangles += be2.rectangles
        be2.rectangles = []
        be3.rectangles += be4.rectangles
        be4.rectangles = []
        be1.bv2 = w
        be3.bv2 = z
        conjugate(be1, be3)
        w.inn.remove(be2)
        z.inn.remove(be4)
        v.inn, v.out = [], []
        y.inn, y.out = [], []
        z.inn.append(be3)
        w.inn.append(be1)
        self.vs.remove(v)
        self.vs.remove(y)
        self.es.remove(be2)
        self.es.remove(be4)

    def condense(self):
        l = len(self.vs)
        for bv in list(self.vs): # copy because can be: "Set changed size during iteration"
            if len(bv.inn) == 1 and len(bv.out) == 1 and (bv.inn[0] != bv.out[0]):
                self.join_biedges(bv.inn[0], bv.out[0])
        print >>logfile, "Condensed %d bi-vertices (out of %d). %d bi-vertices left." % (
            l - len(self.vs), l, len(self.vs))

    def project(self):
        g = Graph()
        for be in self.es:
            # v ---------be--------> w
            # y <-----be.conj------- x
            v, w = be.bv1, be.bv2
            x, y = be.conj.bv1, be.conj.bv2
            g.add_vertex(v.bvid, y.bvid)
            g.add_vertex(w.bvid, x.bvid)
            g.add_vertex(y.bvid, v.bvid)
            g.add_vertex(x.bvid, w.bvid)
            seq = be.get_seq(self.graph.K)
            g.add_edge(be.beid, v.bvid, w.bvid, len(seq) - self.graph.K, be.conj.beid)
            g.add_seq(be.beid, seq)
        g.update_K()
        return g

    def check(self):
        for rectangle in self.rectangles.itervalues():
            assert rectangle.conj, "Some rectangle have no conjugate"
            assert len(rectangle.diagonals) >= 1
            for diagonal in rectangle.diagonals:
                assert diagonal.e1 == rectangle.e1
                assert diagonal.e2 == rectangle.e2
                assert diagonal.e1.len >= diagonal.offsetc
                assert diagonal.e2.len >= diagonal.offsetd
                assert 0 <= diagonal.offseta
                assert 0 <= diagonal.offsetb
                assert diagonal.offsetd - diagonal.offsetb == diagonal.offsetc - diagonal.offseta
        for bvertex in self.vs:
            assert bvertex.conj, "Some bi-vertex have no conjugate"
            assert len(bvertex.inn) == len(bvertex.conj.out)
            assert len(bvertex.out) == len(bvertex.conj.inn)
            for bedge in bvertex.inn:
                assert bedge.conj in bvertex.conj.out
            for bedge in bvertex.out:
                assert bedge.conj in bvertex.conj.inn
        for bedge in self.es:
            assert bedge.conj, "Some bi-edge have no conjugate"

    def filter_diagonals(self):
        before = sum(len(be.rectangles[0].diagonals) for be in self.es)
        for bv in self.vs:
            if len(bv.inn) == 0 or len(bv.out) == 0: # end of contig
                continue
            # check
            for be in itertools.chain(bv.inn, bv.out):
                assert len(be.rectangles) == 1, "Should be run before condensation"
            # get possible in-out keys
            innkeys = set()
            for be in bv.inn:
                innkeys |= be.rectangles[0].keys2
            outkeys = set()
            for be in bv.out:
                outkeys |= be.rectangles[0].keys1
            assert (innkeys & outkeys), "In-rectangles have no shared diagonal with out-rectangles"
            # filter by them
            for be in bv.inn:
                r = be.rectangles[0]
                r.diagonals = filter(lambda d: d.key2 in outkeys, r.diagonals)
                r.update_keys()
            for be in bv.out:
                r = be.rectangles[0]
                r.diagonals = filter(lambda d: d.key1 in innkeys, r.diagonals)
                r.update_keys()
            # check
            for be in itertools.chain(bv.inn, bv.out):
                assert len(be.rectangles[0].diagonals) > 0, "Some rectangle have no diagonal now"
        after = sum(len(be.rectangles[0].diagonals) for be in self.es)
        print >>logfile, "Filtered %d diagonals (out of %d). %d diagonals left." % (before - after, before, after)

########
# MAIN #
########

if __name__ == '__main__':

    ##########
    # PARAMS #
    ##########

    if len(sys.argv) < 5:
        print >>sys.stderr, 'Usage: python', sys.argv[0], 'graph.grp sequences.sqn pairedinfo.prd d [outprefix = out]'
        exit()

    grp_filename = sys.argv[1]
    sqn_filename = sys.argv[2]
    prd_filename = sys.argv[3]
    d = int(sys.argv[4]) # distance between paired reads ( == gap + read_length == insert_size - read_length)
    outprefix = sys.argv[5] if len(sys.argv) > 5 else 'out'

    max_scaffold = d / 2 # TODO: think why?

    logfile = open(outprefix + '.log', 'w') # or sys.stderr

    #########
    # PARSE #
    #########

    ingraph = Graph()
    ingraph.load(grp_filename, sqn_filename)
    ingraph.check()

    #########
    # BUILD #
    #########

    bgraph = BGraph(ingraph)
    bgraph.load(prd_filename)
    bgraph.check()
    bgraph.build()
    bgraph.check()
    bgraph.glue()
    bgraph.check()
    #bgraph.filter_diagonals() # optional
    bgraph.check()
    #bgraph.scaffold() # optional
    ingraph.check()
    bgraph.check()
    bgraph.condense()
    bgraph.check()

    ##########
    # OUTPUT #
    ##########

    outgraph = bgraph.project()
    outgraph.check()
    outgraph.fasta(open(outprefix + '.fasta', 'w'))
    outgraph.stats(logfile)
    outgraph.save(outprefix)