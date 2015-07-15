#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


###
# Note: you should update ../../../cap/cap_environment_manager.hpp before compiling it!
'''
diff --git a/assembler/src/cap/cap_environment_manager.hpp b/assembler/src/cap/cap_environment_manager.hpp
index b452aa7..f68d6c2 100644
--- a/assembler/src/cap/cap_environment_manager.hpp
+++ b/assembler/src/cap/cap_environment_manager.hpp
@@ -4,6 +4,7 @@
 
 #include "compare_standard.hpp"
 #include "graphio.hpp"
+#include "contig_output.hpp"
 
 #include "comparison_utils.hpp"
 #include "diff_masking.hpp"
@@ -284,6 +285,7 @@ class CapEnvironmentManager {
        printer.SaveGraph(filename);
        printer.SaveEdgeSequences(filename);
        printer.SavePositions(filename, *env_->edge_pos_);
+        debruijn_graph::OutputContigs(*env_->graph_, filename + ".fasta");
 
     // Saving coloring of graph
     cap::SaveColoring(*env_->graph_, *env_->int_ids_, *env_->coloring_, filename);
'''
###

import sys
import os
import shutil

sys.path.append(os.path.join(os.path.abspath(sys.path[0]), '../..'))
chain_template = os.path.join(os.path.abspath(sys.path[0]), 'ideal_assembler_chain')
ideal_assembler_bin_dir = os.path.join(os.path.abspath(sys.path[0]), '../../../..')
ideal_assembler_bin = os.path.join(ideal_assembler_bin_dir, 'run_cap')

import fastaparser

def update_template_params(filename, params_subst_dict):
    old_file = open(filename, 'r')
    old_lines = old_file.readlines()
    old_file.close()
    new_file = open(filename, 'w')
    for line in old_lines:
        for k, v in params_subst_dict.items():
            if k in line:
                line = line.replace(k, v)
        new_file.write(line)
    new_file.close()

# MAIN
if len(sys.argv) != 4:
    print("Usage: " + sys.argv[0] + " <input fasta> <K or K1,K2,K3> <output_dir>")
    sys.exit()

if len(sys.argv[2].split(',')) > 1:
    K_list = map(int, sys.argv[2].split(','))
else:
    K_list = [int(sys.argv[2])]
output_dir = os.path.abspath(sys.argv[3])
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

# creating single-entry references and chains
params_subst_dict = dict()
input_fasta = fastaparser.read_fasta(sys.argv[1])
cwd = os.getcwd()
os.chdir(ideal_assembler_bin_dir)
for K in K_list:
    print("Starting with K=" + str(K))
    result_fasta = []
    for id, fasta_entry in enumerate(input_fasta):
        cur_ref_name = os.path.join(output_dir, 'chr_' + str(id) + '.fasta')
        cur_chain_name = os.path.join(output_dir, 'chr_' + str(id) + '_K' + str(K) + '_chain')
        log_filename = os.path.join(output_dir, 'chr_' + str(id) + '_K' + str(K) + '.log')
        fastaparser.write_fasta_to_file(cur_ref_name, [fasta_entry])
        shutil.copy(chain_template, cur_chain_name)
        cur_params_subst_dict = dict(params_subst_dict)
        cur_params_subst_dict['OUT_BASE'] = 'chr_' + str(id) + '_K' + str(K)
        tmp_dir = os.path.join(ideal_assembler_bin_dir, 'data/cap/cache/env_' + cur_params_subst_dict['OUT_BASE'])
        cur_params_subst_dict['REFERENCE'] = cur_ref_name
        cur_params_subst_dict['KMER_SIZE'] = str(K)
        update_template_params(cur_chain_name, cur_params_subst_dict)
        cmd_line = ideal_assembler_bin + ' ' + cur_chain_name + ' >> ' + log_filename + ' 2>> ' +  log_filename
        print('running with ' + os.path.basename(cur_ref_name) + ' on K=' + str(K))
        return_code = os.system(cmd_line)
        if return_code:
            print("Error happened when executing cmd_line " + cmd_line)
            sys.exit(1)
        catch_phrase = 'Outputting contigs to'
        for line in open(log_filename):
            if catch_phrase in line:
                result = os.path.abspath(line.split(catch_phrase)[1].strip())
                break
        print('result is ' + result)
        result_fasta += fastaparser.read_fasta(result)

        os.remove(cur_ref_name)
        os.remove(cur_chain_name)
        os.remove(log_filename)
        shutil.rmtree(tmp_dir)

    # final result
    ideal_assembly_filename = os.path.join(output_dir, os.path.splitext(os.path.basename(sys.argv[1]))[0] + '_K' + str(K) + '.fasta')
    print("Final result for K = " + str(K) + " is " + ideal_assembly_filename)
    fastaparser.write_fasta_to_file(ideal_assembly_filename, result_fasta)
print("Finished")

    


