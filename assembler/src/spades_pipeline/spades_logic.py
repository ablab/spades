#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import shutil
import glob

import support
import process_cfg
from process_cfg import bool_to_str
from process_cfg import load_config_from_file

def prepare_config_spades(filename, cfg, log, prev_K, K, last_one):
    subst_dict = dict()

    subst_dict["K"] = str(K)
    subst_dict["run_mode"] = "false"
    subst_dict["dataset"] = cfg.dataset
    subst_dict["output_base"] = cfg.output_dir
    subst_dict["additional_contigs"] = cfg.additional_contigs
    subst_dict["entry_point"] = "construction"
    subst_dict["developer_mode"] = bool_to_str(cfg.developer_mode)
    subst_dict["SAM_writer_enable"] = bool_to_str(cfg.generate_sam_files and last_one)
    subst_dict["align_original_reads"] = bool_to_str(cfg.align_original_reads)
    subst_dict["align_before_RR"] = bool_to_str(not cfg.paired_mode)
    subst_dict["align_after_RR"] = bool_to_str(cfg.paired_mode)
    subst_dict["gap_closer_enable"] = bool_to_str(last_one and cfg.gap_closer)
    subst_dict["paired_mode"] = bool_to_str(last_one and cfg.paired_mode)
    subst_dict["additional_ec_removing"] = bool_to_str(last_one)
    subst_dict["use_additional_contigs"] = bool_to_str(prev_K)
    subst_dict["max_threads"] = cfg.max_threads
    subst_dict["max_memory"] = cfg.max_memory
    subst_dict["correct_mismatches"] = bool_to_str(last_one)
    if "resolving_mode" in cfg.__dict__:
        subst_dict["resolving_mode"] = cfg.resolving_mode

    process_cfg.substitute_params(filename, subst_dict, log)

def get_read_length(output_dir, K):
    estimated_params = load_config_from_file(os.path.join(output_dir, "K%d" % (K), "_est_params.info"))
    return estimated_params.__dict__["RL"]

def run_iteration(configs_dir, execution_home, cfg, log, K, use_additional_contigs, last_one):
    data_dir = os.path.join(cfg.output_dir, "K%d" % (K))
    if os.path.exists(data_dir):
        shutil.rmtree(data_dir)
    os.makedirs(data_dir)
    bin_reads_dir = os.path.join(cfg.output_dir, ".bin_reads")

    dst_configs = os.path.join(data_dir, "configs")
    shutil.copytree(os.path.join(configs_dir, "debruijn"), dst_configs)
    cfg_file_name = os.path.join(dst_configs, "config.info")
    # removing template configs
    for root, dirs, files in os.walk(dst_configs):
        for cfg_file in files:
            cfg_file = os.path.join(root, cfg_file)
            if cfg_file.endswith('.info.template'):
                if os.path.isfile(cfg_file.split('.template')[0]):
                    os.remove(cfg_file)
                else:
                    os.rename(cfg_file, cfg_file.split('.template')[0])

    prepare_config_spades(cfg_file_name, cfg, log, use_additional_contigs, K, last_one)
    prev_K = K

    command = os.path.join(execution_home, "spades") + " " +\
               os.path.abspath(cfg_file_name)

    if os.path.isdir(bin_reads_dir):
        if glob.glob(os.path.join(bin_reads_dir, "*_cor*")):
            for cor_filename in glob.glob(os.path.join(bin_reads_dir, "*_cor*")):
                cor_index = cor_filename.rfind("_cor")
                new_bin_filename = cor_filename[:cor_index] + cor_filename[cor_index + 4:]
                shutil.move(cor_filename, new_bin_filename)

    log.info("\n== Running assembler: " + ("K%d" % (K)) + "\n")

    support.sys_call(command, log, execution_home)

def run_spades(configs_dir, execution_home, cfg, log):
    if not isinstance(cfg.iterative_K, list):
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    bin_reads_dir = os.path.join(cfg.output_dir, ".bin_reads")
    if os.path.isdir(bin_reads_dir):
        shutil.rmtree(bin_reads_dir)

    if len(cfg.iterative_K) == 1:
        run_iteration(configs_dir, execution_home, cfg, log, cfg.iterative_K[0], False, True)
        K = cfg.iterative_K[0]
    else:
        run_iteration(configs_dir, execution_home, cfg, log, cfg.iterative_K[0], False, False)
        RL = get_read_length(cfg.output_dir, cfg.iterative_K[0])
        if (cfg.iterative_K[1] > RL):
            if cfg.paired_mode:
                log.info("Second value of iterative K exceeded estimated read length. Rerunning in paired mode for the first value of K")
                run_iteration(configs_dir, execution_home, cfg, log, cfg.iterative_K[0], False, True)
                K = cfg.iterative_K[0]
        else:
            rest_of_iterative_K = cfg.iterative_K
            rest_of_iterative_K.pop(0)
            count = 0
            for K in rest_of_iterative_K:
                count += 1
                last_one = count == len(cfg.iterative_K) or rest_of_iterative_K[count] > RL
                run_iteration(configs_dir, execution_home, cfg, log, K, True, last_one)
                if last_one:
                    break
            if count < len(cfg.iterative_K):
                log.info("Iterations stopped. Value of K exceeded estimated read length")

    latest = os.path.join(cfg.output_dir, "K%d" % (K))

    if os.path.isfile(os.path.join(latest, "final_contigs.fasta")):    
        shutil.copyfile(os.path.join(latest, "final_contigs.fasta"), cfg.result_contigs)
    if cfg.paired_mode:
        if os.path.isfile(os.path.join(latest, "scaffolds.fasta")):
            shutil.copyfile(os.path.join(latest, "scaffolds.fasta"), cfg.result_scaffolds)

    if cfg.developer_mode:
        # before repeat resolver contigs
        # before_RR_contigs = os.path.join(os.path.dirname(cfg.result_contigs), "simplified_contigs.fasta")
        # shutil.copyfile(os.path.join(latest, "simplified_contigs.fasta"), before_RR_contigs)
        # saves
        saves_link = os.path.join(os.path.dirname(cfg.result_contigs), "saves")
        if os.path.lexists(saves_link): # exists return False for broken link! lexists return True
            os.remove(saves_link)
        os.symlink(os.path.join(latest, "saves"), saves_link)

    #    os.remove(cfg.additional_contigs)

    #if glob.glob(os.path.join(latest, "*.sam")):
    #    sam_file_linkname = os.path.join(os.path.dirname(cfg.result_contigs),
    #        "contigs.sam")
    #    if os.path.exists(sam_file_linkname):
    #        os.remove(sam_file_linkname)
    #    os.symlink(glob.glob(os.path.join(latest, "*.sam"))[0], sam_file_linkname)

    if os.path.isdir(bin_reads_dir):
        shutil.rmtree(bin_reads_dir)
    if not cfg.paired_mode:
        return cfg.result_contigs, "", latest

    return cfg.result_contigs, cfg.result_scaffolds, latest
