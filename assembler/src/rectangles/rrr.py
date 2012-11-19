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
import os
import utils
import fastaparser
import saveparser
from test_util import TestUtils
from graph import Graph
from bigraph import BGraph
from rectangle_set import RectangleSet
import experimental
import logging
import check_diags
from optparse import OptionParser

def makelogger(logfilename):
    # create logger with 'rectangles'
    logger = logging.getLogger('rectangles')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(logfilename, mode='w')
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)


def resolve(input_path, output_path, test_utils, genome, is_sc):

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    grp_filename = os.path.join(input_path, 'late_pair_info_counted.grp')
    sqn_filename = os.path.join(input_path, 'late_pair_info_counted.sqn')
    cvr_filename = os.path.join(input_path, 'late_pair_info_counted.cvr')
    first_prd_filename = os.path.join(input_path, 'late_pair_info_counted.prd')
    prd_filename = os.path.join(input_path, 'late_pair_info_counted.prd' if experimental.filter != experimental.Filter.spades else 'distance_filling_cl.prd')
    pst_filename = os.path.join(input_path, 'distance_estimation.pst') if experimental.filter == experimental.Filter.pathsets else None
    inf_filename = os.path.join(input_path, 'late_pair_info_counted_est_params.info')
    log_filename = os.path.join(output_path, 'rectangles.log')
    config = saveparser.config(inf_filename)
    d = config.median - config.RL

    makelogger(log_filename)
    logger = logging.getLogger('rectangles')

    logger.info("Rectangle Resolving %s..." % input_path)
    logger.info("d = %d..." % d)

    #################################
    # PARSE INITIAL BE BRUIJN GRAPH #
    #################################

    ingraph = Graph()
    ingraph.load(grp_filename, sqn_filename, cvr_filename)
    ingraph.check()
    edges_before_loop_DG = ingraph.find_loops(10, 1000) 
    maxN50 = 0
    maxgraph = None
    maxbgraph = None
    maxthreshold = 0

    rs = RectangleSet(ingraph, d, test_utils, prd_filename, first_prd_filename, config)
    if experimental.filter == experimental.Filter.pathsets:
        rs.pathsets(pst_filename)
    else:
        rs.filter(prd_filename, config)
    logger.info("  RectangleSet built.")
    if experimental.filter == experimental.Filter.spades:
        thresholds = [0.0] # everything supported by paired info
    elif experimental.filter == experimental.Filter.pathsets:
        thresholds = [-1] # everything from pathsets file
    else:
        thresholds = rs.percentiles()
    logger.info("  Checking thresholds %s..." % thresholds)
    for threshold in set(thresholds):
        logger.info("  Checking threshold %f..." % threshold)
        bgraph = rs.bgraph(threshold)
        if not bgraph.diagonals:
            continue
        bgraph.build_missing_rectangles(ingraph.K, rs)
        bgraph.condense()
        outgraph = bgraph.project(output_path, is_sc)
        maxgraph = outgraph
        maxbgraph = bgraph
        maxthreshold = threshold

    maxgraph.fasta(open(os.path.join(output_path, 'begin_rectangles.fasta'), 'w'))
    #maxgraph.save(os.path.join(output_path, 'rectangles'))
    maxbgraph.save(output_path, ingraph.K)
    maxbgraph.check_tips(ingraph.K)
    outgraph = maxbgraph.project(output_path, is_sc)
    outgraph.fasta(open(os.path.join(output_path, 'delete_tips.fasta'), 'w'))
    edges_before_loop = maxbgraph.delete_loops(ingraph.K, 1000, 10)
    maxbgraph.condense()
    outgraph = maxbgraph.project(output_path, is_sc)
    outgraph.fasta(open(os.path.join(output_path,"delete_tips_delete_loops_1000.fasta"),"w"))
    to_del = set()
    for eid in edges_before_loop_DG:
          if eid in edges_before_loop:
            to_del.add(eid)
    print "to_del", len(to_del)
    for eid in to_del:
      del edges_before_loop_DG[eid]
    maxbgraph.delete_missing_loops(edges_before_loop_DG, ingraph.K, 1000, 10)
    maxbgraph.condense()
    outgraph = maxbgraph.project(output_path, is_sc)
    outgraph.fasta(open(os.path.join(output_path, 'delete_tips_delete_all_loops_1000.fasta'), 'w'))
    edges_before_loop_DG = ingraph.find_loops(4, 10000) 
    edges_before_loop_DG = edges_before_loop_DG or maxbgraph.delete_missing_loops(ingraph.K, 10000,10)
    to_del = set()
    for eid in edges_before_loop_DG:
          if eid in edges_before_loop:
            to_del.add(eid)
    print "to_del", len(to_del)
    for eid in to_del:
      del edges_before_loop_DG[eid]
    maxbgraph.delete_missing_loops(edges_before_loop_DG, ingraph.K, 10000, 10)
    maxbgraph.condense()
    outgraph = maxbgraph.project(output_path, is_sc)
    outgraph.fasta(open(os.path.join(output_path, "after_deleting_big_loops.fasta"), "w"))
    additional_paired_info = dict()
    should_connect = maxbgraph.edges_expand(5000)
    should_connect_by_first_pair_info = maxbgraph.use_scaffold_paired_info(2 * maxbgraph.d, rs.additional_prd)
    for (e1id, e2id) in should_connect_by_first_pair_info:
      if e1id not in additional_paired_info and maxbgraph.es[e1id].conj.eid not in additional_paired_info and e2id not in additional_paired_info:
        additional_paired_info[e1id] = [maxbgraph.es[e1id], maxbgraph.es[e2id]]
        additional_paired_info[maxbgraph.es[e1id].conj.eid] = [maxbgraph.es[e2id].conj, maxbgraph.es[e1id].conj]
    outgraph.fasta_for_long_contigs(ingraph.K, maxbgraph.d, is_sc, open(os.path.join(output_path,"rectangles_extend.fasta"),"w"), should_connect, additional_paired_info)
    outgraph.fasta_for_long_contigs(ingraph.K, maxbgraph.d, is_sc, open(os.path.join(output_path,"rectangles_extend_before_scaffold.fasta"),"w"), should_connect, dict())
    
    maxbgraph.print_about_edges([20586, 23014, 23806, 19630,23350], ingraph.K)
    outgraph.save(os.path.join(output_path,"last_graph"))
    if genome:  
      check_diags.check(genome, maxbgraph, maxgraph.K, open(os.path.join(output_path, "check_log.txt"), "w"), test_utils) 

    logger.info("Best Threshold = %d" % maxthreshold)
    logger.info("Best N50 = %d" % maxN50)

def parser_options():
    parser = OptionParser()
    parser.add_option("-s", "", dest="saves_dir", help="Name of directory with saves")
    parser.add_option("-o", "", dest="out_dir", help = "Output directory, default = out (optional)", default = "out")
    parser.add_option("-g", "", dest ="genome", help = "File with genome (optional)") 
    parser.add_option("-d", "", dest = "debug_logger", help = "File for debug logger (optional)", default = "debug_log.txt")
    parser.add_option("-k", "", type = int, dest = "k", help = "k (optional)")
    parser.add_option("-D", "", type = int, dest="d", help = "d (optional)")
    parser.add_option("", "--sc", dest = "sc", action="store_true", help = "Turn on if data is sincle-cell (optional)", default = False)
    return parser

def make_rectangles_from_genome(options):
    k = options.k
    ingraph = Graph()
    _, genome = fastaparser.read_fasta(options.genome).next()
    ingraph.make_graph(genome, int(k))
    edges_before_loop_DG = ingraph.find_loops(10, 1000) 
    ingraph.save(os.path.join(options.out_dir,"graph"))
    rs = RectangleSet(ingraph, int(options.d))
    rs.filter_without_prd()
    f_left = open(os.path.join(options.out_dir, "paired_genom_contigs_1.fasta"),"w") # TODO: what is it?
    f_right = open(os.path.join(options.out_dir, "paired_genom_contigs_2.fasta"),"w") # TODO: what is it?
    contigs_id = 0
    for key, rect in rs.rectangles.items():
      for key, diag in rect.diagonals.items():
        e1 = rect.e1.seq
        e2 = rect.e2.seq
        f_left.write(">" + str(contigs_id) + "/1\n")
        f_left.write(e1[diag.offseta:diag.offsetc])
        f_left.write("\n")
        f_right.write(">"+str(contigs_id) + "/2\n")
        f_right.write(e2[diag.offsetb:diag.offsetd])
        f_right.write("\n")
        contigs_id += 1
    bgraph = rs.bgraph_from_genome()
    bgraph.condense()
    outgraph = bgraph.project(options.out_dir, False)
    outgraph.fasta(open(os.path.join(options.out_dir, 'rectangles.fasta'), 'w'))
    
if __name__ == '__main__':

    ##########
    # PARAMS #
    ##########
    parser = parser_options()
    (options, args) = parser.parse_args()

    if not os.path.exists(options.out_dir):
        os.mkdir(options.out_dir)
    
    if options.genome and not options.saves_dir:  
      if not options.k or not options.d:
        print "specify k and d"
        sys.exit(1)
      make_rectangles_from_genome(options)
      sys.exit(1)

    if not options.saves_dir:
      parser.print_help()
      sys.exit(0)
    
    input_dir = options.saves_dir
    outpath = options.out_dir
    reference_information_file = os.path.join(input_dir,"late_pair_info_counted_etalon_distance.txt")
    test_util = TestUtils(reference_information_file, os.path.join(outpath, options.debug_logger))
    resolve(input_dir, outpath, test_util, options.genome, options.sc)
    
    if test_util.has_ref_info:
      test_util.stat()
