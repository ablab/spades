#!/usr/bin/env python

import os
import shutil
import glob

import support
import process_cfg
from process_cfg import bool_to_str

def prepare_config_spades(filename, cfg, prev_K, K, last_one):
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

    process_cfg.substitute_params(filename, subst_dict)


def run_spades(spades_home, execution_home, cfg):
    if not isinstance(cfg.iterative_K, list):
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    count = 0
    prev_K = None

    bin_reads_dir = os.path.join(cfg.output_dir, ".bin_reads")
    if os.path.isdir(bin_reads_dir):
        shutil.rmtree(bin_reads_dir)

    for K in cfg.iterative_K:
        count += 1

        dst_configs = os.path.join(cfg.output_dir, "K%d" % (K))
        if os.path.exists(dst_configs):
            shutil.rmtree(dst_configs)
        os.makedirs(dst_configs)

        dst_configs = os.path.join(dst_configs, "configs")
        shutil.copytree(os.path.join(spades_home, "configs", "debruijn"), dst_configs)
        cfg_file_name = os.path.join(dst_configs, "config.info")
        # removing template configs
        for root, dirs, files in os.walk(dst_configs):
            for cfg_file in files:
                if cfg_file.endswith('.template'):
                    os.remove(os.path.join(root, cfg_file))

        prepare_config_spades(cfg_file_name, cfg, prev_K, K,
            count == len(cfg.iterative_K))
        prev_K = K

        command = ""
        if "use_jemalloc" in cfg.__dict__ and os.path.isfile("jemalloc.sh"):
            command = os.path.abspath("jemalloc.sh") + " "

        command += os.path.join(execution_home, "spades") + " " +\
                   os.path.abspath(cfg_file_name)

        if os.path.isdir(bin_reads_dir):
            if glob.glob(os.path.join(bin_reads_dir, "*_cor*")):
                for cor_filename in glob.glob(os.path.join(bin_reads_dir, "*_cor*")):
                    cor_index = cor_filename.rfind("_cor")
                    new_bin_filename = cor_filename[:cor_index] + cor_filename[cor_index + 4:]
                    shutil.move(cor_filename, new_bin_filename)

        print("\n== Running assembler: " + command + "\n")

        support.sys_call(command, execution_home)

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
        if os.path.exists(saves_link):
            os.remove(saves_link)
        os.symlink(os.path.join(latest, "saves"), saves_link)

    #    os.remove(cfg.additional_contigs)

    if glob.glob(os.path.join(latest, "*.sam")):
    #        sam_file_linkname = os.path.join(os.path.dirname(cfg.result_contigs),
    #            "contigs.sam")
        if os.path.exists(sam_file_linkname):
            os.remove(sam_file_linkname)
        #        os.symlink(glob.glob(os.path.join(latest, "*.sam"))[0], sam_file_linkname)

    if os.path.isdir(bin_reads_dir):
        shutil.rmtree(bin_reads_dir)
    if not cfg.paired_mode:
        return cfg.result_contigs, "", latest

    return cfg.result_contigs, cfg.result_scaffolds, latest
