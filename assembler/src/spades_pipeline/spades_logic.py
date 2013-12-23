#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import shutil
import support
import process_cfg
from process_cfg import bool_to_str
from site import addsitedir
import options_storage

BASE_STAGE = "construction"
READS_TYPES_USED_IN_CONSTRUCTION = ["paired-end", "single"]

def prepare_config_spades(filename, cfg, log, additional_contigs_fname, K, stage, saves_dir, last_one):
    subst_dict = dict()

    subst_dict["K"] = str(K)
    subst_dict["run_mode"] = "false"
    if "diploid_mode" in cfg.__dict__:
        subst_dict["diploid_mode"] = bool_to_str(cfg.diploid_mode)
    subst_dict["dataset"] = process_cfg.process_spaces(cfg.dataset)
    subst_dict["output_base"] = process_cfg.process_spaces(cfg.output_dir)
    subst_dict["tmp_dir"] = process_cfg.process_spaces(cfg.tmp_dir)
    if additional_contigs_fname:
        subst_dict["additional_contigs"] = process_cfg.process_spaces(additional_contigs_fname)
        subst_dict["use_additional_contigs"] = bool_to_str(True)
    else:
        subst_dict["use_additional_contigs"] = bool_to_str(False)
    subst_dict["main_iteration"] = bool_to_str(last_one)
    subst_dict["entry_point"] = stage
    subst_dict["load_from"] = saves_dir
    subst_dict["developer_mode"] = bool_to_str(cfg.developer_mode)
    subst_dict["gap_closer_enable"] = bool_to_str(last_one)
    subst_dict["rr_enable"] = bool_to_str(last_one and cfg.rr_enable)
#    subst_dict["topology_simplif_enabled"] = bool_to_str(last_one)
    subst_dict["max_threads"] = cfg.max_threads
    subst_dict["max_memory"] = cfg.max_memory
    subst_dict["correct_mismatches"] = bool_to_str(last_one)
    if "resolving_mode" in cfg.__dict__:
        subst_dict["resolving_mode"] = cfg.resolving_mode
    if "careful" in cfg.__dict__:
        subst_dict["mismatch_careful"] = bool_to_str(cfg.careful)
    if "pacbio_mode" in cfg.__dict__:
        subst_dict["pacbio_test_on"] = bool_to_str(cfg.pacbio_mode)
        subst_dict["pacbio_reads"] = process_cfg.process_spaces(cfg.pacbio_reads)

    process_cfg.substitute_params(filename, subst_dict, log)


def get_read_length(output_dir, K, ext_python_modules_home, log):
    est_params_filename = os.path.join(output_dir, "K%d" % K, "final.lib_data")
    max_read_length = 0
    if os.path.isfile(est_params_filename):
        addsitedir(ext_python_modules_home)
        if sys.version.startswith('2.'):
            import pyyaml2 as pyyaml
        elif sys.version.startswith('3.'):
            import pyyaml3 as pyyaml
        est_params_data = pyyaml.load(open(est_params_filename, 'r'))
        for reads_library in est_params_data:
            if reads_library['type'] in READS_TYPES_USED_IN_CONSTRUCTION:
                if int(reads_library["read length"]) > max_read_length:
                    max_read_length = int(reads_library["read length"])
    if max_read_length == 0:
        support.error("Failed to estimate maximum read length! File with estimated params: " + est_params_filename, log)
    return max_read_length


def update_k_mers_in_special_cases(cur_k_mers, RL, log, silent=False):
    if options_storage.auto_K_allowed():
        if RL >= 250:
            if not silent:
                support.warning("Default k-mer sizes were set to %s because estimated "
                                "read length (%d) is equal to or greater than 250" % (str(options_storage.K_MERS_250), RL), log)
            return options_storage.K_MERS_250
        if RL >= 150:
            if not silent:
                support.warning("Default k-mer sizes were set to %s because estimated "
                                "read length (%d) is equal to or greater than 150" % (str(options_storage.K_MERS_150), RL), log)
            return options_storage.K_MERS_150
    return cur_k_mers


def run_iteration(configs_dir, execution_home, cfg, log, K, prev_K, last_one):
    data_dir = os.path.join(cfg.output_dir, "K%d" % K)
    stage = BASE_STAGE
    saves_dir = os.path.join(data_dir, 'saves')
    dst_configs = os.path.join(data_dir, "configs")
    cfg_file_name = os.path.join(dst_configs, "config.info")

    if options_storage.continue_mode:
        if os.path.isfile(os.path.join(data_dir, "final_contigs.fasta")) and not (options_storage.restart_from and
            (options_storage.restart_from == ("k%d" % K) or options_storage.restart_from.startswith("k%d:" % K))):
            log.info("\n== Skipping assembler: " + ("K%d" % K) + " (already processed)")
            return
        if options_storage.restart_from and options_storage.restart_from.find(":") != -1:
            stage = options_storage.restart_from[options_storage.restart_from.find(":") + 1:]
        support.continue_from_here(log)

    if stage != BASE_STAGE:
        if not os.path.isdir(saves_dir):
            support.error("Cannot restart from stage %s: saves were not found (%s)!" % (stage, saves_dir))
    else:
        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)
        os.makedirs(data_dir)

        shutil.copytree(os.path.join(configs_dir, "debruijn"), dst_configs)
        # removing template configs
        for root, dirs, files in os.walk(dst_configs):
            for cfg_file in files:
                cfg_file = os.path.join(root, cfg_file)
                if cfg_file.endswith('.info.template'):
                    if os.path.isfile(cfg_file.split('.template')[0]):
                        os.remove(cfg_file)
                    else:
                        os.rename(cfg_file, cfg_file.split('.template')[0])

    log.info("\n== Running assembler: " + ("K%d" % K) + "\n")
    if prev_K:
        additional_contigs_fname = os.path.join(cfg.output_dir, "K%d" % prev_K, "simplified_contigs.fasta")
        if not os.path.isfile(additional_contigs_fname):
            support.warning("additional contigs for K=%d were not found (%s)!" % (K, additional_contigs_fname), log)
            additional_contigs_fname = None
    else:
        additional_contigs_fname = None
    if "read_buffer_size" in cfg.__dict__:
        construction_cfg_file_name = os.path.join(dst_configs, "construction.info")
        process_cfg.substitute_params(construction_cfg_file_name, {"read_buffer_size": cfg.read_buffer_size}, log)
    prepare_config_spades(cfg_file_name, cfg, log, additional_contigs_fname, K, stage, saves_dir, last_one)

    command = [os.path.join(execution_home, "spades"), cfg_file_name]

## this code makes sense for src/debruijn/simplification.cpp: corrected_and_save_reads() function which is not used now
#    bin_reads_dir = os.path.join(cfg.output_dir, ".bin_reads")
#    if os.path.isdir(bin_reads_dir):
#        if glob.glob(os.path.join(bin_reads_dir, "*_cor*")):
#            for cor_filename in glob.glob(os.path.join(bin_reads_dir, "*_cor*")):
#                cor_index = cor_filename.rfind("_cor")
#                new_bin_filename = cor_filename[:cor_index] + cor_filename[cor_index + 4:]
#                shutil.move(cor_filename, new_bin_filename)
    support.sys_call(command, log)


def run_spades(configs_dir, execution_home, cfg, dataset_data, ext_python_modules_home, log):
    if not isinstance(cfg.iterative_K, list):
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    # checking and removing conflicting K-mer directories
    if options_storage.restart_from:
        processed_K = []
        for k in range(options_storage.MIN_K, options_storage.MAX_K, 2):
            cur_K_dir = os.path.join(cfg.output_dir, "K%d" % k)
            if os.path.isdir(cur_K_dir) and os.path.isfile(os.path.join(cur_K_dir, "final_contigs.fasta")):
                processed_K.append(k)
        if processed_K:
            RL = get_read_length(cfg.output_dir, processed_K[0], ext_python_modules_home, log)
            needed_K = update_k_mers_in_special_cases(cfg.iterative_K, RL, log, silent=True)
            needed_K = [k for k in needed_K if k < RL]
            k_to_delete = []
            for id, k in enumerate(needed_K):
                if len(processed_K) == id:
                    k_to_delete = [processed_K[-1]] # the last K in processed K was run in "last_one" mode
                    break
                if processed_K[id] != k:
                    k_to_delete = processed_K[id:]
                    break
            if not k_to_delete and (len(processed_K) > len(needed_K)):
                k_to_delete = processed_K[len(needed_K) - 1:]
            if k_to_delete:
                log.info("Restart mode: removing previously processed directories for K=%s "
                         "to except conflicts with K specified with --restart-from" % (str(k_to_delete)))
                for k in k_to_delete:
                    shutil.rmtree(os.path.join(cfg.output_dir, "K%d" % k))

    bin_reads_dir = os.path.join(cfg.output_dir, ".bin_reads")
    if os.path.isdir(bin_reads_dir) and not options_storage.continue_mode:
        shutil.rmtree(bin_reads_dir)
    cfg.tmp_dir = support.get_tmp_dir(prefix="spades_")

    if len(cfg.iterative_K) == 1:
        run_iteration(configs_dir, execution_home, cfg, log, cfg.iterative_K[0], None, True)
        K = cfg.iterative_K[0]
    else:
        run_iteration(configs_dir, execution_home, cfg, log, cfg.iterative_K[0], None, False)
        prev_K = cfg.iterative_K[0]
        RL = get_read_length(cfg.output_dir, cfg.iterative_K[0], ext_python_modules_home, log)
        cfg.iterative_K = update_k_mers_in_special_cases(cfg.iterative_K, RL, log)
        if cfg.iterative_K[1] + 1 > RL:
            if cfg.rr_enable:
                support.warning("Second value of iterative K (%d) exceeded estimated read length (%d). "
                                "Rerunning for the first value of K (%d) with Repeat Resolving" %
                                (cfg.iterative_K[1], RL, cfg.iterative_K[0]), log)
                run_iteration(configs_dir, execution_home, cfg, log, cfg.iterative_K[0], None, True)
                K = cfg.iterative_K[0]
        else:
            rest_of_iterative_K = cfg.iterative_K
            rest_of_iterative_K.pop(0)
            count = 0
            for K in rest_of_iterative_K:
                count += 1
                last_one = count == len(cfg.iterative_K) or (rest_of_iterative_K[count] + 1 > RL)
                run_iteration(configs_dir, execution_home, cfg, log, K, prev_K, last_one)
                prev_K = K
                if last_one:
                    break
            if count < len(cfg.iterative_K):
                support.warning("Iterations stopped. Value of K (%d) exceeded estimated read length (%d)" %
                                (cfg.iterative_K[count], RL), log)

    latest = os.path.join(cfg.output_dir, "K%d" % K)

    for format in [".fasta", ".fastg"]:
        if os.path.isfile(os.path.join(latest, "before_rr" + format)):
            result_before_rr_contigs = os.path.join(os.path.dirname(cfg.result_contigs), "before_rr" + format)
            if not os.path.isfile(result_before_rr_contigs) or not options_storage.continue_mode:
                shutil.copyfile(os.path.join(latest, "before_rr" + format), result_before_rr_contigs)
        if os.path.isfile(os.path.join(latest, "final_contigs" + format)):
            if not os.path.isfile(cfg.result_contigs[:-6] + format) or not options_storage.continue_mode:
                shutil.copyfile(os.path.join(latest, "final_contigs" + format), cfg.result_contigs[:-6] + format)
        if cfg.rr_enable:
            if os.path.isfile(os.path.join(latest, "scaffolds" + format)):
                if not os.path.isfile(cfg.result_scaffolds[:-6] + format) or not options_storage.continue_mode:
                    shutil.copyfile(os.path.join(latest, "scaffolds" + format), cfg.result_scaffolds[:-6] + format)

    if cfg.developer_mode:
        # saves
        saves_link = os.path.join(os.path.dirname(cfg.result_contigs), "saves")
        if os.path.lexists(saves_link): # exists return False for broken link! lexists return True
            os.remove(saves_link)
        os.symlink(os.path.join(latest, "saves"), saves_link)

    if os.path.isdir(bin_reads_dir):
        shutil.rmtree(bin_reads_dir)
    if os.path.isdir(cfg.tmp_dir):
        shutil.rmtree(cfg.tmp_dir)

    return latest
