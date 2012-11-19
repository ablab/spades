#######################################################
###
### Parsing module for .grp / .sqn / .prd files
###
#######################################################

from collections import namedtuple
import utils

########################
# PARSE GRAPH VERTICES #
########################

def grp_vertices(filename):
    grp = open(filename)
    V, E = map(int, grp.readline().split())
    for _ in xrange(V): # Vertex 2 ~ 1 .
        vertex, vid, tilda, compl, dot = grp.readline().split()
        assert vertex == 'Vertex', "grp format violation"
        assert tilda == '~', "grp format violation"
        assert dot == '.', "grp format violation"
        vid = int(vid)
        compl = int(compl)
        yield vid, compl
    grp.close()

#####################
# PARSE GRAPH EDGES #
#####################

# parse edges from grp file
# yields: (eid, v1id, v2id, l, compl)
def grp_edges(filename):
    grp = open(filename)
    V, E = map(int, grp.readline().split())
    for _ in xrange(V): # Vertex 2 ~ 1 .
        grp.readline()
    line = grp.readline() # empty line
    assert line.strip() == '', "grp format violation"
    for _ in xrange(E): # Edge 15 : 12 -> 2, l = 838 ~ 16 .
        edge, eid, colon, v1id, arrow, v2id, letterl, equal, l, tilda, compl, dot = grp.readline().split()
        assert edge == 'Edge', "grp format violation"
        assert colon == ':', "grp format violation"
        assert arrow == '->', "grp format violation"
        assert v2id[-1] == ',', "grp format violation"
        assert letterl == 'l', "grp format violation"
        assert equal == '=', "grp format violation"
        assert tilda == '~', "grp format violation"
        assert dot == '.', "grp format violation"
        eid = int(eid)
        v1id = int(v1id)
        v2id = int(v2id[:-1])
        l = int(l)
        compl = int(compl)
        yield eid, v1id, v2id, l, compl
    grp.close()

#####################
# PARSE PAIRED INFO #
#####################

# parse paired info from prd file
# yields: (e1id, e2id, D, weight, delta)
def prd(filename):
    prd = open(filename)
    P = int(prd.readline())
    for _ in xrange(P): # 15 15 0.00 1.91 0.00 .
        e1id, e2id, dist, weight, delta, dot = prd.readline().split()
        assert dot == '.', "prd format violation"
        D = float(dist)
        if D < 0:
            continue
        e1id = int(e1id)
        e2id = int(e2id)
        weight = float(weight)
        delta = float(delta)
        yield e1id, e2id, D, weight, delta
    prd.close()

##################
# PARSE PATHSETS #
##################

# parse paired info from prd file
# yields: e1, e2, [ [path1], [path2], ..., [pathN]
def pst(filename):
    pst = open(filename)
    P = int(pst.readline()) # 154
    for _ in xrange(P): # 1  834
        X = int(pst.readline().split()[0])
        pathset = []
        e1 = None
        e2 = None
        for _ in xrange(X): # 3  :  2659 4516 7237
            line = pst.readline().split()
            size = int(line[0])
            assert line[1] == ':', "pst format violation"
            assert len(line) == size + 2, "pst format violation"
            path = [int(x) for x in line[2:]]
            assert not e1 or e1 == path[0]
            assert not e2 or e2 == path[-1]
            e1 = path[0]
            e2 = path[-1]
            pathset.append(path)
        if X: # vpyl to work with Dima's IS480
            assert e1 and e2, (e1, e2)
            yield e1, e2, pathset
        pst.readline()
    pst.close()

###################
# PARSE SEQUENCES #
###################

# parse sequences from sqn file
# yields: (eid, seq)
def sqn(filename):
    sqn = open(filename)
    while True: #15 AGCTTTTCATTCTG .
        eid = sqn.readline()
        seq = sqn.readline().strip()
        if (not eid):
            break
        assert eid[0] == '>', "sqn format violation (it's just FASTA!)"
        eid = int(eid[1:])
        yield eid, seq
    sqn.close()


##################
# PARSE COVERAGE #
##################

# parse coverage from cvr file
# yields: (eid, cvr)
def cvr(filename):
    cvrf = open(filename)
    E = int(cvrf.readline()) # 154
    for line in cvrf:
        eid, cvr, dot = line.split()
        assert dot == '.'
        eid = int(eid)
        cvr = float(cvr)
        yield eid, cvr
    cvrf.close()

#######################
# PARSE CONFIG (info) #
#######################

Config = namedtuple('Config', ['RL', 'IS', 'is_var', 'perc', 'avg_coverage', 'median', 'mad', 'hist', 'K'])

def config(filename):
    f = open(filename)
    RL = int(f.readline().split()[1])
    IS = float(f.readline().split()[1])
    is_var = float(f.readline().split()[1])
    perc = ' '.join(f.readline().split()[1:])
    avg_coverage = float(f.readline().split()[1])
    median = int(f.readline().split()[1])
    mad = int(f.readline().split()[1])
    hist = ' '.join(f.readline().split()[1:])
    try:
        K = int(f.readline().split()[1])
    except: # for old saves
        K = 55
    f.close()
    hist = utils.strtomap(hist)
    hist2 = {}
    for k, v in hist.iteritems():
        hist2[k - RL] = v # * 1.0 / s
    return Config(RL, IS, is_var, perc, avg_coverage, median, mad, hist2, K)

#######
# EOF #
#######
