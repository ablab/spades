#!/usr/bin/python
#
# calculate fp, fn, tp, tn for diagonals, by reference

import os
import utils

def check(reference, bgraph, K, log, test_util):
    K = 15
    #logstream = open(os.path.join(folder, 'rectangles.log'), 'a')
    #corr = os.path.join(folder, 'rectangles.corr.diagonals')
    #wrong = os.path.join(folder, 'rectangles.wrong.diagonals')

    max_mismatch = 2 # for match()

    # read reference genome 
    ref = open(reference)
    header = ref.readline()
    assert header[0] == '>'
    genome = ''
    for line in ref:
        genome += line.strip()
    ref.close()
    rcgenome = utils.rc(genome)

    def get_index(genome, not_rc):
        s = genome + genome[:K]
        d = {}
        for i in xrange(len(genome)):
            kmer = s[i:i+K]
            if kmer in d:
                d[kmer].add(i+1 if not_rc else -(len(genome)-i - K + 1))
            else:
                d[kmer] = {i+1 if not_rc else -(len(genome) -i - K + 1)}
        return d

    igenome = get_index(genome, True)
    ircgenome = get_index(rcgenome, False)


    def find_index_begin(seq):
       ref_indexs =set()
       for i in range(len(seq) - K + 1):  
          kmer = seq[i:i + K]
          if kmer in igenome:
            ref_index = igenome[kmer]
            for ind in ref_index:
              ref_indexs.add(ind - i)
          if kmer in ircgenome:
            ref_index = ircgenome[kmer]
            for ind in ref_index:
              ref_indexs.add(ind - i)
       return ref_indexs  

    count_true_be = 0
    count_true_etalon_dist_be = 0
    count_part_true_be = 0
    count_part_true_etalon_dist_be = 0
    count_be = 0
    covered = []
    for be in bgraph.bes:
      count_be +=1 
      is_true_be = True
      is_part_true = False
      is_true_etalon_dist = True
      is_part_true_etalon_dist = False
      for diag in be.diagonals:
        is_true_diag = test_util.is_true_diagonal(diag)
        if is_true_diag:
          is_part_true_etalon_dist = True
        else:
          is_true_etalon_dist = False
        seq1 = diag.rectangle.e1.seq[diag.offseta:diag.offsetc + 55]
        seq2 = diag.rectangle.e2.seq[diag.offsetb:diag.offsetd + 55]
        indexs1 = find_index_begin(seq1)
        indexs2 = find_index_begin(seq2)
        is_true = False
        log.write(str( indexs1)+"\n")
        log.write(str( indexs2)+"\n")
        log.write(str( diag.d)+"\n")
        for i1 in indexs1:
            for i2 in indexs2:
              if diag.d-5<=abs(i1 - i2) <= diag.d+5:
                is_true = True
                covered.append((i1, i1+len(seq1))) 
        if is_true:
          is_part_true = True
        else:
          is_true_be = False
        log.write("diag is true " + str(is_true) +"\n")
      if is_part_true and not is_true_be:
        count_part_true_be +=1
      if is_true_be:
        count_true_be +=1
      if is_part_true_etalon_dist and not is_true_etalon_dist: 
        count_part_true_etalon_dist_be += 1
      if is_true_etalon_dist:
        count_true_etalon_dist_be += 1

      #print diag, count_part_true_be, count_true_be
    log.write( "true " + str(count_true_be)+ str( " part_true ") + str( count_part_true_be)+ " false " + str( count_be - count_true_be - count_part_true_be) + "\n")
    log.write( "etalon true " + str(count_true_etalon_dist_be)+ str( " part_true ") + str( count_part_true_etalon_dist_be)+ " false " + str( count_be - count_true_etalon_dist_be - count_part_true_etalon_dist_be) + "\n")
    
    covered.sort()
   # print covered
    missed = 0
    end = covered[0][1]
    for cov in covered:
      if cov[0] > end:
        missed += 1
      end = cov[1]
    log.write("missing" + str( missed) +"\n")
    
    
    
    """def match(genome, seq, pos):
        l = len(genome)
        mism = 0
        for c in seq:
            if genome[pos] != c:
                mism += 1
                if mism > max_mismatch: return False
            pos += 1
            if pos == l: pos = 0
        return True

    def check(genome, igenome, seq1, seq2, d):
        kmer1 = seq1[:K]
        if kmer1 not in igenome: return False
        kmer2 = seq2[:K]
        if kmer2 not in igenome: return False
        set1 = igenome[kmer1]
        set2 = igenome[kmer2]
        for pos in set1:
            if pos+d in set2:
                if match(genome, seq1, pos) and match(genome, seq2, pos+d):
                    return True
        return False

    def align(filename, genome, igenome, rcgenome, ircgenome):
        f = open(filename)
        d = int(f.readline())
        true, false = 0, 0
        for cnt, line in enumerate(f):
            seq1, seq2 = line.split()
            ok = check(genome, igenome, seq1, seq2, d) or check(rcgenome, ircgenome, seq1, seq2, d)
            if ok: true += 1
            else: false += 1
            if cnt % 100 == 0: print >>logstream, filename, cnt, (true, false)
        f.close()
        return true, false


    tp, fp = align(corr, genome, igenome, rcgenome, ircgenome)
    fn, tn = align(wrong, genome, igenome, rcgenome, ircgenome)

    # http://en.wikipedia.org/wiki/Sensitivity_and_specificity
    sensitivity = tp * 1.0 / (tp + fn)
    specificity = tn * 1.0 / (fp + tn)
    positive_predicted_value = tp * 1.0 / (tp + fp)
    negative_predicted_value = tn * 1.0 / (fn + tn)

    precision = positive_predicted_value
    recall = sensitivity
    F = 2.0 * precision * recall / (precision + recall)

    print >>logstream, '-----'
    print >>logstream, 'TP = %5d | FP = %5d' % (tp, fp)
    print >>logstream, 'FN = %5d | TN = %5d' % (fn, tn)
    print >>logstream, '-----'
    print >>logstream, 'Sensitivity = %f' % sensitivity
    print >>logstream, 'Specificity = %f' % specificity
    print >>logstream, '-----'
    print >>logstream, 'Positive predictive value = %f' % positive_predicted_value
    print >>logstream, 'Negative predictive value = %f' % negative_predicted_value
    print >>logstream, '-----'
    print >>logstream, 'Recall (same as PPV) = %f' % recall
    print >>logstream, 'Precision (same as sensitivity) = %f' % precision
    print >>logstream, '-----'
    print >>logstream, 'F-measure = %f' % F"""
