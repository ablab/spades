
#######################################################
###
### Parsing module for .grp / .sqn / .prd files
###
#######################################################

########################
# PARSE GRAPH VERTICES #
########################

# parse vertices from grp file
# yields: (vid, compl)
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
# yields: (e1id, e2id, D, weight, var)
def prd(filename):
    prd = open(filename)
    P = int(prd.readline())
    for _ in xrange(P): # 15 15 0.00 1.91 0.00 .
        e1id, e2id, dist, weight, var, dot = prd.readline().split()
        assert dot == '.', "prd format violation"
        D = float(dist)
        assert D*2 == int(round(D*2)), "prd format violation: distance should be integer or half-integer"
        if D < 0:
            continue
        e1id = int(e1id)
        e2id = int(e2id)
        weight = float(weight)
        var = float(var)
        Ds = range(int(round(D-var)), int(round(D+var)) + 1)
        yield e1id, e2id, Ds, weight
    prd.close()


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

###################
# PARSE SEQUENCES #
###################

def rc(seq):
    return seq.translate('*****************************************************************TVGHEFCDIJMLKNOPQYSAUBWXRZ[\]^_`tvghefcdijmlknopqysaubwxrz*************************************************************************************************************************************')[::-1]

#######
# EOF #
#######
