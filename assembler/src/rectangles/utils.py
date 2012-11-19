##################
# UTIL FUNCTIONS #
##################

#import math

def conjugate(a, b):
    a.conj = b
    if b: b.conj = a

def seq_join(seq1, seq2): # join two sequences, prioritize ACGT over N
    if len(seq1) != len(seq2):
      print seq1
      print seq2
    assert len(seq1) == len(seq2)
#    for pos, (c1, c2) in enumerate(itertools.izip(seq1, seq2)):
#        assert c1 == c2 or c1 == 'N' or c2 == 'N', (c1, c2, pos, len(seq1), seq1[pos-20:pos+20], seq2[pos-20:pos+20]) # seq1 + '|' + seq2,
#    seq = ''.join(map(lambda x, y: x if x != 'N' else (y if y != 'N' else 'N'), seq1, seq2))
#    return seq
    seq1 = seq1.strip('N')
    return seq1
    d = len(seq2) - len(seq1)
    seq = ""
    different = 0
    for i in range(len(seq1)):
      if seq1[i] == 'N':
        seq += seq2[i]
      else:
        seq += seq1[i]
      if seq1[i] != 'N' and seq2[i] != 'N' and seq1[i]!= seq2[i]:
        different += 1
    if different > 0:
      print "different\n",seq1, "\n", seq2
    #elif len(seq1) > 4*d:
    # seq += seq2[len(seq1):]
    #elif len(seq1) < d:
    #  return seq
    #else:
    #  return seq[d:]

    return seq

def erf(x): # from http://stackoverflow.com/a/457805/92396
    return 0
    #return math.erf(x)

def rc(seq):
    return seq.translate('*****************************************************************TVGHEFCDIJMLKNOPQYSAUBWXRZ[\]^_`tvghefcdijmlknopqysaubwxrz*************************************************************************************************************************************')[::-1]

def strtomap(str):
    # e.g
    # "195 28650;196 36394;197 46938;198 57458;199 70565"
    # to
    # {195: 28650, 196 : 36394, 197 : 46938, 198 : 57458, 199 : 70565}
    str = str.strip('"')
    res = {}
    for pair in str.split(';'):
        x, y = map(int, pair.split())
        res[x] = y
    return res


def quantile(ls, q):
    assert 0 <= q <= 100
    q *= 0.01
    pos = int(round((len(ls)-1) * q))
    return ls[pos]

def choose_d(config):
    #return int(round(config.IS)) - config.RL
    mx = (-1, None)
    interval = 2*(config.RL - config.K) # to maximize sum of histogram over this interval
    for k, v in config.hist.iteritems():
        s = 0
        for x in xrange(interval):
            if (k+x) in config.hist:
                s += config.hist[k+x]
        mx = max(mx, (s, k))
    d = mx[1] + interval / 2
    return d

def get_overlap(seq1, seq2):
  overlap = 0
  while len(seq1) - 1 - overlap >= 0 and overlap< len(seq2) and seq1[len(seq1) - 1 - overlap] == seq2[overlap]:
    overlap += 1
  return overlap
