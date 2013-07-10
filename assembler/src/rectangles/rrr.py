#!/usr/bin/python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


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


def save_fasta(bgraph, output_path, is_sc, name_file, is_careful):
    bgraph.condense()
    graph = bgraph.project(output_path, is_sc, is_careful)
    graph.fasta(open(os.path.join(output_path, name_file), 'w'))
    return graph

def resolve(input_path, output_path, test_utils, genome, is_sc, is_careful):

    grp_filename = os.path.join(input_path, 'late_pair_info_counted.grp')
    sqn_filename = os.path.join(input_path, 'late_pair_info_counted.sqn')
    cvr_filename = os.path.join(input_path, 'late_pair_info_counted.cvr')
    first_prd_filename = os.path.join(input_path, 'late_pair_info_counted_0.prd')

    if experimental.filter != experimental.Filter.spades:
        prd_filename = first_prd_filename
    else:
        prd_filename = os.path.join(input_path, 'distance_estimation_0_cl.prd')
    if experimental.filter == experimental.Filter.pathsets:
        pst_filename = os.path.join(input_path, 'distance_estimation.pst')
    inf_filename = os.path.join(input_path, 'late_pair_info_counted_est_params.info')
    log_filename = os.path.join(output_path, 'rectangles.log')
    config = saveparser.config(inf_filename)

    d = config['median'] - config['RL']
    
    if d <= 0:
      print "Read length", config ['RL'], "is smaller than insert size", config['median'],", can't do anything"
      return
    
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
    logger.info("init rectangles set")
    rs = RectangleSet(ingraph, d, test_utils, prd_filename, first_prd_filename, config)
    if experimental.filter == experimental.Filter.pathsets:
        rs.pathsets(pst_filename)
    else:
        logger.info("begin filter")
        rs.filter(prd_filename, config)
    logger.info("  RectangleSet built.")
    
    threshold = 0.0
    logger.info("  Checking threshold %f..." % threshold)
    maxbgraph = rs.bgraph(threshold)
    save_fasta(maxbgraph, output_path, is_sc, 'begin_rectangles.fasta', is_careful)
    logger.info("outputed begin rectangles")
    maxbgraph.check_tips(ingraph.K)
    save_fasta(maxbgraph, output_path, is_sc, 'delete_tips.fasta', is_careful)
    logger.info("outputed delete tips")
    edges_before_loop = maxbgraph.delete_loops(ingraph.K, 1000, 10)
    save_fasta(maxbgraph, output_path, is_sc, "delete_tips_delete_loops_1000.fasta", is_careful)
    logger.info("outputed delete loops")
    edges_before_loop_DG = ingraph.find_loops(10, 1000, rs)
    logger.info("find DG 1000 loops")
    to_del = set(edges_before_loop_DG.keys()) & edges_before_loop
    for eid in to_del:
        del edges_before_loop_DG[eid]
    
    maxbgraph.delete_missing_loops(edges_before_loop_DG, ingraph.K, 1000, 10)
    logger.info("delete missing loops")
    save_fasta(maxbgraph, output_path, is_sc, 'delete_tips_delete_all_loops_1000.fasta', is_careful)
    
    edges_before_loop_DG = ingraph.find_loops(4, 10000, rs)
    to_del = set(edges_before_loop_DG.keys()) & edges_before_loop
    for eid in to_del:
        del edges_before_loop_DG[eid]
    
    maxbgraph.delete_missing_loops(edges_before_loop_DG, ingraph.K, 10000, 10)
    outgraph = save_fasta(maxbgraph, output_path, is_sc, "after_deleting_big_loops.fasta", is_careful)

    additional_paired_info = dict()
    should_connect = maxbgraph.edges_expand(5000)
    should_connect_by_first_pair_info = maxbgraph.use_scaffold_paired_info(2 * maxbgraph.d, rs.prd_for_scaffold)
    
    for (e1id, e2id) in should_connect_by_first_pair_info:
        if e1id not in additional_paired_info and maxbgraph.es[
                                                  e1id].conj.eid not in additional_paired_info and e2id not in additional_paired_info:
            additional_paired_info[e1id] = [maxbgraph.es[e1id], maxbgraph.es[e2id]]
            additional_paired_info[maxbgraph.es[e1id].conj.eid] = [maxbgraph.es[e2id].conj, maxbgraph.es[e1id].conj]
    
    outgraph.fasta_for_long_contigs(ingraph.K, maxbgraph.d, is_sc, is_careful,
        open(os.path.join(output_path, "rectangles_extend.fasta"), "w"), should_connect, additional_paired_info)
    outgraph.fasta_for_long_contigs(ingraph.K, maxbgraph.d, is_sc, is_careful,
        open(os.path.join(output_path, "rectangles_extend_before_scaffold.fasta"), "w"), should_connect, dict())

    outgraph.save(os.path.join(output_path, "last_graph"))
    
    if genome:
        check_diags.check(genome, maxbgraph, maxgraph.K, open(os.path.join(output_path, "check_log.txt"), "w"),
            test_utils)

def parser_options():
    parser = OptionParser()
    parser.add_option("-s", "", dest="saves_dir", help="Name of directory with saves")
    parser.add_option("-g", "", dest="genome", help="File with genome (optional)")
    parser.add_option("-o", "", dest="out_dir", help="Output directory, default = out (optional)", default="out")
    parser.add_option("-d", "", dest="debug_logger", help="File for debug logger (optional)", default="debug_log.txt")
    parser.add_option("", "--sc", dest="sc", action="store_true", help="Turn on if data is sincle-cell (optional)",
        default=False)
    parser.add_option("", "--careful", dest="is_careful", action="store_true", help="make contigs careful, can reduce N50 and genes", default= False)
    return parser


if __name__ == '__main__':
    ##########
    # PARAMS #
    ##########
    parser = parser_options()
    (options, args) = parser.parse_args()

    if not os.path.exists(options.out_dir):
        os.mkdir(options.out_dir)

    if not options.saves_dir:
        parser.print_help()
        sys.exit(0)

    input_dir = options.saves_dir
    outpath = options.out_dir
    reference_information_file = os.path.join(input_dir, "late_pair_info_counted_etalon_distance.txt")
    test_util = TestUtils(reference_information_file, os.path.join(outpath, options.debug_logger))
    resolve(input_dir, outpath, test_util, options.genome, options.sc, options.is_careful)

    if test_util.has_ref_info:
        test_util.stat()
