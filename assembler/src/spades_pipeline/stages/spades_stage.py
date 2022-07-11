#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
import sys
from site import addsitedir

import commands_parser
import options_storage
from stages import stage
from stages import spades_iteration_stage
import support
from process_cfg import merge_configs


def get_read_length(output_dir, K, ext_python_modules_home, log):
    est_params_filename = os.path.join(output_dir, "K%d" % K, "final.lib_data")
    max_read_length = 0
    if os.path.isfile(est_params_filename):
        addsitedir(ext_python_modules_home)
        if sys.version.startswith("2."):
            import pyyaml2 as pyyaml
        elif sys.version.startswith("3."):
            import pyyaml3 as pyyaml
        est_params_data = pyyaml.load(open(est_params_filename))
        max_read_length = int(est_params_data["nomerge max read length"])
        log.info("Max read length detected as %d" % max_read_length)
    if max_read_length == 0:
        support.error("Failed to estimate maximum read length! File with estimated params: %s" % est_params_filename, log)
    return max_read_length


def update_k_mers_in_special_cases(cur_k_mers, RL, log, silent=False):
    if options_storage.auto_K_allowed():
        if RL >= 250:
            if not silent:
                log.info("Default k-mer sizes were set to %s because estimated "
                         "read length (%d) is equal to or greater than 250" % (str(options_storage.K_MERS_250), RL))
            return options_storage.K_MERS_250
        if RL >= 150:
            if not silent:
                log.info("Default k-mer sizes were set to %s because estimated "
                         "read length (%d) is equal to or greater than 150" % (str(options_storage.K_MERS_150), RL))
            return options_storage.K_MERS_150
    if RL <= max(cur_k_mers):
        new_k_mers = [k for k in cur_k_mers if k < RL]
        if not silent:
            log.info("K-mer sizes were set to %s because estimated "
                     "read length (%d) is less than %d" % (str(new_k_mers), RL, max(cur_k_mers)))
        return new_k_mers
    return cur_k_mers


def reveal_original_k_mers(RL):
    if options_storage.original_k_mers is None or options_storage.original_k_mers == "auto":
        cur_k_mers = options_storage.args.k_mers
        options_storage.args.k_mers = options_storage.original_k_mers
        original_k_mers = update_k_mers_in_special_cases(options_storage.K_MERS_SHORT, RL, None, silent=True)
        options_storage.args.k_mers = cur_k_mers
    else:
        original_k_mers = options_storage.original_k_mers
    original_k_mers = [k for k in original_k_mers if k < RL]
    return original_k_mers


def rna_k_values(support, dataset_data, log):
    rna_rl = support.get_reads_length(dataset_data, log, ["merged reads"])
    upper_k = int(rna_rl / 2) - 1
    if upper_k % 2 == 0:
        upper_k -= 1

    lower_k = min(max(int(rna_rl / 3), options_storage.RNA_MIN_K), options_storage.RNA_MAX_LOWER_K)
    if lower_k % 2 == 0:
        lower_k -= 1

    use_iterative = True
    if upper_k <= lower_k:
        use_iterative = False

    if upper_k < options_storage.RNA_MIN_K:
        support.warning("\nauto K value (%d) is too small, recommended to be at least %d.\n" % (upper_k, options_storage.RNA_MIN_K))
        if rna_rl <= options_storage.RNA_MIN_K:
            support.warning(
                "read length is too small (%d), but keeping current K value anyway. Consider setting K manually. K\n" % (
                    rna_rl))
        else:
            upper_k = options_storage.RNA_MIN_K
        log.info("Upper K value is set to %d.\n" % (upper_k))

    if upper_k > options_storage.MAX_K:
        log.info("\nAuto K value (%d) is too large, all K values should not exceed %d. Setting k=%d.\n"
                 % (upper_k, options_storage.MAX_K, options_storage.MAX_K))
        upper_k = options_storage.MAX_K

    if not use_iterative:
        return [upper_k]
    return [lower_k, upper_k]


def generateK_for_rna(cfg, dataset_data, log):
    if cfg.iterative_K == "auto":
        k_values = options_storage.K_MERS_RNA
        if not options_storage.args.iontorrent:
            k_values = rna_k_values(support, dataset_data, log)
        cfg.iterative_K = k_values
        log.info("K values to be used: " + str(k_values))


def generateK_for_rnaviral(cfg, dataset_data, log):
    if cfg.iterative_K == "auto":
        k_values = rna_k_values(support, dataset_data, log)
        # FIXME: Hack-hack-hack! :)
        if min(k_values) == options_storage.RNA_MAX_LOWER_K:
            k_values = [options_storage.K_MERS_RNA[0]] + k_values
        if min(k_values) > 21:
            k_values = [21] + k_values
        cfg.iterative_K = k_values
        log.info("K values to be used: " + str(k_values))


def generateK(cfg, log, dataset_data, silent=False):
    if options_storage.args.rna:
        generateK_for_rna(cfg, dataset_data, log)
    elif options_storage.args.rnaviral:
        generateK_for_rnaviral(cfg, dataset_data, log)
    elif not options_storage.args.iontorrent:
        RL = support.get_primary_max_reads_length(dataset_data, log, ["merged reads"],
                                                  options_storage.READS_TYPES_USED_IN_CONSTRUCTION)
        if options_storage.auto_K_allowed():
            if options_storage.args.plasmid:
                if RL >= 150:
                    if not silent:
                        log.info("Default k-mer sizes were set to %s because estimated read length (%d) is equal to or greater than 150" % (str(options_storage.K_MERS_PLASMID_LONG), RL))
                    cfg.iterative_K = options_storage.K_MERS_PLASMID_LONG
                else:
                    if not silent:
                        log.info("Default k-mer sizes were set to %s because estimated read length (%d) is less than 150" % (str(options_storage.K_MERS_PLASMID_100), RL))
                    cfg.iterative_K = options_storage.K_MERS_PLASMID_100
            else:
                if RL >= 250:
                    if not silent:
                        log.info("Default k-mer sizes were set to %s because estimated "
                                 "read length (%d) is equal to or greater than 250" % (str(options_storage.K_MERS_250), RL))
                    cfg.iterative_K = options_storage.K_MERS_250
                elif RL >= 150:
                    if not silent:
                        log.info("Default k-mer sizes were set to %s because estimated "
                                 "read length (%d) is equal to or greater than 150" % (str(options_storage.K_MERS_150), RL))
                    cfg.iterative_K = options_storage.K_MERS_150
        if RL <= max(cfg.iterative_K):
            new_k_mers = [k for k in cfg.iterative_K if k < RL]
            if not silent:
                log.info("K-mer sizes were set to %s because estimated "
                         "read length (%d) is less than %d" % (str(new_k_mers), RL, max(cfg.iterative_K)))
            cfg.iterative_K = new_k_mers

    if not isinstance(cfg.iterative_K, list):
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

class PlasmidGlueFileStage(stage.Stage):
    STAGE_NAME = "metaextrachromosomal glue files"
    def __init__(self, latest, *args):
        super(PlasmidGlueFileStage, self).__init__(*args)
        self.latest = latest

    def get_command(self, cfg):
        self.cfg = cfg
        args = [os.path.join(self.python_modules_home, "spades_pipeline", "scripts", "plasmid_glue.py")]
        args.append(self.latest)
        command = [commands_parser.Command(STAGE=self.STAGE_NAME,
                                       path=sys.executable,
                                       args=args,
                                       short_name=self.short_name,
                                       )]
        return command 


class SpadesCopyFileStage(stage.Stage):
    STAGE_NAME = "Copy files"

    def always_copy(self, output_file, latest, cfg):
        return not cfg.correct_scaffolds

    def rna_copy(self, output_file, latest, cfg):
        return options_storage.args.rna and self.always_copy(output_file, latest, cfg)

    def has_hmm(self, output_file = None, latest = None, cfg = None):
        return options_storage.args.bio or options_storage.args.custom_hmms or options_storage.args.corona

    def correct_scaffolds_copy(self, output_file, latest, cfg):
        return cfg.correct_scaffolds

    def not_rna_copy(self, output_file, latest, cfg):
        return (not options_storage.args.rna) and self.always_copy(output_file, latest, cfg)

    def rr_enably_copy(self, output_file, latest, cfg):
        return (not options_storage.args.rna) and cfg.rr_enable and self.always_copy(output_file, latest, cfg)

    class OutputFile(object):
        def __init__(self, output_file, tmp_file, need_to_copy):
            self.output_file = output_file
            self.tmp_file = tmp_file
            self.need_to_copy = need_to_copy

    def set_output_files(self):

        self.output = [
            self.OutputFile(self.cfg.result_scaffolds, "corrected_scaffolds.fasta", self.correct_scaffolds_copy),
            self.OutputFile(os.path.join(os.path.dirname(self.cfg.result_contigs), "before_rr.fasta"),
                            "before_rr.fasta", self.always_copy),
            self.OutputFile(os.path.join(os.path.dirname(self.cfg.result_contigs), "assembly_graph_after_simplification.gfa"),
                            "assembly_graph_after_simplification.gfa", self.always_copy),
            self.OutputFile(self.cfg.result_transcripts, "transcripts.fasta", self.rna_copy),
            self.OutputFile(self.cfg.result_transcripts_paths, "transcripts.paths", self.rna_copy),
            self.OutputFile(self.cfg.result_contigs, "final_contigs.fasta", self.not_rna_copy),
            self.OutputFile(os.path.join(os.path.dirname(self.cfg.result_contigs),
                                         "first_pe_contigs.fasta"), "first_pe_contigs.fasta", self.not_rna_copy),
            self.OutputFile(os.path.join(os.path.dirname(self.cfg.result_contigs),
                                         "strain_graph.gfa"), "strain_graph.gfa", self.not_rna_copy),
            self.OutputFile(self.cfg.result_scaffolds, "scaffolds.fasta", self.rr_enably_copy),
            self.OutputFile(self.cfg.result_scaffolds_paths, "scaffolds.paths", self.rr_enably_copy),
            self.OutputFile(self.cfg.result_graph_gfa, "assembly_graph_with_scaffolds.gfa", self.always_copy),
            self.OutputFile(self.cfg.result_graph, "assembly_graph.fastg", self.always_copy),
            self.OutputFile(self.cfg.result_contigs_paths, "final_contigs.paths", self.not_rna_copy),
            self.OutputFile(self.cfg.result_gene_clusters, "gene_clusters.fasta", self.has_hmm),
            self.OutputFile(self.cfg.result_bgc_statistics, "hmm_statistics.txt", self.has_hmm),
            self.OutputFile(self.cfg.result_domain_graph, "domain_graph.dot", self.has_hmm)
        ]

        for filtering_type in options_storage.filtering_types:
            prefix = filtering_type + "_filtered_"
            result_filtered_transcripts = os.path.join(self.cfg.output_dir,
                                                       prefix + options_storage.transcripts_name)
            self.output.append(self.OutputFile(result_filtered_transcripts, prefix + "final_paths.fasta", self.rna_copy))

    def __init__(self, latest, *args):
        super(SpadesCopyFileStage, self).__init__(*args)
        self.latest = latest

    def copy_files(self):
        latest = self.latest
        for outputfile in self.output:
            if outputfile.need_to_copy(outputfile, latest, self.cfg):
                filename = os.path.join(latest, outputfile.tmp_file)
                if os.path.isfile(filename):
                    shutil.copyfile(filename, outputfile.output_file)

    def get_command(self, cfg):
        self.cfg = cfg
        self.set_output_files()
        args = [os.path.join(self.python_modules_home, "spades_pipeline", "scripts", "copy_files.py")]
        for outputfile in self.output:
            if outputfile.need_to_copy(outputfile, self.latest, self.cfg):
                filename = os.path.join(self.latest, outputfile.tmp_file)
                args.append(filename)
                args.append(outputfile.output_file)
        bin_reads_dir = os.path.join(self.cfg.output_dir, ".bin_reads")
        command = [commands_parser.Command(STAGE=self.STAGE_NAME,
                                           path=sys.executable,
                                           args=args,
                                           short_name=self.short_name,
                                           del_after=[os.path.relpath(bin_reads_dir, options_storage.args.output_dir),
                                                      os.path.relpath(self.cfg.tmp_dir, options_storage.args.output_dir)])]
        return command


class SpadesStage(stage.Stage):
    stages = []
    STAGE_NAME = "Assembling"
    latest = ""

    def __init__(self, cfg, get_stage, *args):
        super(SpadesStage, self).__init__(*args)
        self.get_stage = get_stage

        self.generate_cfg(cfg, self.output_files)

        # creating dataset
        dataset_filename = os.path.join(self.cfg.output_dir, "dataset.info")
        with open(dataset_filename, 'w') as dataset_file:
            import process_cfg
            # TODO don't exists at that moment, what better to write????
            if self.output_files["corrected_dataset_yaml_filename"] != "":
                dataset_file.write("reads\t%s\n" %
                                   process_cfg.process_spaces(self.output_files["corrected_dataset_yaml_filename"]))
            else:
                dataset_file.write(
                    "reads\t%s\n" % process_cfg.process_spaces(cfg["dataset"].yaml_filename))
            if self.cfg.developer_mode and "reference" in cfg["dataset"].__dict__:
                dataset_file.write("reference_genome\t")
                dataset_file.write(process_cfg.process_spaces(cfg["dataset"].reference) + '\n')

        if not os.path.isdir(self.output_files["misc_dir"]):
            os.makedirs(self.output_files["misc_dir"])

        generateK(self.cfg, self.log, self.dataset_data)

        self.used_K = []
        count = 0
        prev_K = None
        for K in self.cfg.iterative_K:
            count += 1
            last_one = count == len(self.cfg.iterative_K)

            iter_stage = spades_iteration_stage.IterationStage(K, prev_K, last_one, self.get_stage, self.latest,
                                                               "k%d" % K,
                                                               self.output_files, self.tmp_configs_dir,
                                                               self.dataset_data, self.log, self.bin_home,
                                                               self.ext_python_modules_home,
                                                               self.python_modules_home)
            self.stages.append(iter_stage)
            self.latest = os.path.join(self.cfg.output_dir, "K%d" % K)

            self.used_K.append(K)
            prev_K = K
            if last_one:
                break

        if options_storage.args.plasmid and options_storage.args.meta:
            self.stages.append(PlasmidGlueFileStage(self.latest, "plasmid_copy_files", 
                                                    self.output_files,
                                                    self.tmp_configs_dir,
                                                    self.dataset_data, self.log,
                                                    self.bin_home,
                                                    self.ext_python_modules_home,
                                                    self.python_modules_home))
        self.stages.append(SpadesCopyFileStage(self.latest, "copy_files",
                                               self.output_files,
                                               self.tmp_configs_dir,
                                               self.dataset_data, self.log,
                                               self.bin_home,
                                               self.ext_python_modules_home,
                                               self.python_modules_home))

    def generate_cfg(self, cfg, output_files):
        self.cfg = merge_configs(cfg["assembly"], cfg["common"])
        self.cfg.__dict__["result_contigs"] = output_files["result_contigs_filename"]
        self.cfg.__dict__["result_scaffolds"] = output_files["result_scaffolds_filename"]
        self.cfg.__dict__["result_graph"] = output_files["result_assembly_graph_filename"]
        self.cfg.__dict__["result_graph_gfa"] = output_files["result_assembly_graph_filename_gfa"]
        self.cfg.__dict__["result_contigs_paths"] = output_files["result_contigs_paths_filename"]
        self.cfg.__dict__["result_scaffolds_paths"] = output_files["result_scaffolds_paths_filename"]
        self.cfg.__dict__["result_transcripts"] = output_files["result_transcripts_filename"]
        self.cfg.__dict__["result_transcripts_paths"] = output_files["result_transcripts_paths_filename"]
        self.cfg.__dict__["result_gene_clusters"] = output_files["result_gene_clusters_filename"]
        self.cfg.__dict__["result_bgc_statistics"] = output_files["result_bgc_stats_filename"]
        self.cfg.__dict__["result_domain_graph"] = output_files["result_domain_graph_filename"]

        if self.cfg.disable_rr:
            self.cfg.__dict__["rr_enable"] = False
        else:
            self.cfg.__dict__["rr_enable"] = True

        dataset_filename = os.path.join(self.cfg.output_dir, "dataset.info")
        self.cfg.__dict__["dataset"] = dataset_filename
        self.cfg.tmp_dir = support.get_tmp_dir(prefix="spades_")

    def generate_config(self, cfg):
        for stage in self.stages:
            stage.generate_config(self.cfg)

    def get_command(self, cfg):
        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name + "_start")] + \
               [x for stage in self.stages for x in stage.get_command(self.cfg)] + \
               [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name + "_finish")]


def add_to_pipeline(pipeline, get_stage, cfg, output_files, tmp_configs_dir, dataset_data, log, bin_home,
                    ext_python_modules_home, python_modules_home):
    if "assembly" in cfg:
        pipeline.add(SpadesStage(cfg, get_stage, "as", output_files, tmp_configs_dir,
                                 dataset_data, log, bin_home, ext_python_modules_home, python_modules_home))
