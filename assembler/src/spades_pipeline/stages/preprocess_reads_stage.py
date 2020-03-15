#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import gzip

import support
import options_storage
import commands_parser
from stages import stage


class PreprocessInterlacedReads(stage.Stage):
    STAGE_NAME = "Preprocess interlaced reads"

    # {input_filename, out_left_filename, out_right_filename, was_compressed, is_fastq}}
    update_list = []

    def split_interlaced_reads(self, dataset_data, dst, log):
        self.dst = dst
        for reads_library in dataset_data:
            copy_reads_library = dict(reads_library)
            for key, value in copy_reads_library.items():
                if key == "interlaced reads":
                    if "left reads" not in reads_library:
                        reads_library["left reads"] = []
                        reads_library["right reads"] = []
                    for interlaced_reads in value:
                        if interlaced_reads in options_storage.dict_of_prefixes:
                            ext = options_storage.dict_of_prefixes[interlaced_reads]
                        else:
                            ext = os.path.splitext(interlaced_reads)[1]
                        was_compressed = False
                        if ext.endswith(".gz"):
                            was_compressed = True
                            ungzipped = os.path.splitext(interlaced_reads)[0]
                            out_basename, ext = os.path.splitext(os.path.basename(ungzipped))
                        else:
                            out_basename, ext = os.path.splitext(os.path.basename(interlaced_reads))

                        if interlaced_reads in options_storage.dict_of_prefixes:
                            ext = options_storage.dict_of_prefixes[interlaced_reads]

                        if ext.lower().startswith(".fq") or ext.lower().startswith(".fastq"):
                            is_fastq = True
                            ext = ".fastq"
                        else:
                            is_fastq = False
                            ext = ".fasta"

                        out_left_filename = os.path.join(dst, "%s_1%s" % (out_basename, ext))
                        out_right_filename = os.path.join(dst, "%s_2%s" % (out_basename, ext))

                        self.update_list.append({"input_filename": interlaced_reads,
                                                 "out_left_filename": out_left_filename,
                                                 "out_right_filename": out_right_filename,
                                                 "was_compressed": was_compressed,
                                                 "is_fastq": is_fastq})

                        reads_library["left reads"].append(out_left_filename)
                        reads_library["right reads"].append(out_right_filename)

                        if interlaced_reads in options_storage.dict_of_prefixes:
                            del options_storage.dict_of_prefixes[interlaced_reads]
                    del reads_library["interlaced reads"]

    def generate_config(self, cfg):
        self.split_interlaced_reads(self.dataset_data, self.dst, self.log)

        with open(os.path.join(self.tmp_dir, "interlaced"), "w") as fw:
            for update_item in self.update_list:
                fw.write(update_item["input_filename"] + "\n")
                fw.write(update_item["out_left_filename"] + "\n")
                fw.write(update_item["out_right_filename"] + "\n")
                fw.write(str(update_item["was_compressed"]) + "\n")
                fw.write(str(update_item["is_fastq"]) + "\n")

    def get_command(self, cfg):
        command = [commands_parser.Command(STAGE=self.STAGE_NAME,
                                           path=sys.executable,
                                           args=[
                                               os.path.join(self.python_modules_home, "spades_pipeline", "scripts",
                                                            "preprocess_interlaced_reads.py"),
                                               "--args_filename", os.path.join(self.tmp_dir, "interlaced"),
                                               "--dst", self.dst],
                                           short_name=self.short_name)]
        return command

    def __init__(self, dir_for_split_reads, tmp_dir, *args):
        super(PreprocessInterlacedReads, self).__init__(*args)
        self.dst = dir_for_split_reads
        self.tmp_dir = tmp_dir


class PreprocessNxmateReads(stage.Stage):
    STAGE_NAME = "Preprocess nxmate reads"

    # {infile1, infile2}
    update_list = []

    def get_new_names(self, infilename1, infilename2):
        basename1 = os.path.basename(infilename1)
        if os.path.splitext(basename1)[1] == ".gz":
            basename1 = os.path.splitext(basename1)[0]
        basename2 = os.path.basename(infilename2)
        if os.path.splitext(basename2)[1] == ".gz":
            basename2 = os.path.splitext(basename2)[0]
        # open three outfiles
        splitfilenameleft = os.path.join(self.dst, "R1_IJS7_mates_ICC4_" + basename1)

        splitfilenameright = os.path.join(self.dst, "R2_IJS7_mates_ICC4_" + basename2)

        unsplitfilename = os.path.join(self.dst, "unsplit_IJS7_mates_ICC4_" + basename1.replace("_R1_", "_R1R2_"))
        return splitfilenameleft, splitfilenameright, unsplitfilename

    def process_nxmate_reads(self, dataset_data, dst, log):
        self.dst = dst
        for reads_library in dataset_data:
            if reads_library["type"] == "nxmate":
                raw_left_reads = reads_library["left reads"]
                raw_right_reads = reads_library["right reads"]
                reads_library["left reads"] = []
                reads_library["right reads"] = []
                reads_library["single reads"] = []
                for id, left_reads_fpath in enumerate(raw_left_reads):
                    right_reads_fpath = raw_right_reads[id]
                    processed_left_reads_fpath, processed_right_reads_fpath, single_reads_fpath = \
                        self.get_new_names(left_reads_fpath, right_reads_fpath)
                    reads_library["left reads"].append(processed_left_reads_fpath)
                    reads_library["right reads"].append(processed_right_reads_fpath)
                    reads_library["single reads"].append(single_reads_fpath)
                    self.update_list.append({"infile1": left_reads_fpath,
                                             "infile2": right_reads_fpath})
                reads_library["type"] = "mate-pairs"
                reads_library["orientation"] = "fr"

    def generate_config(self, cfg):
        self.process_nxmate_reads(self.dataset_data, self.dst, self.log)
        with open(os.path.join(self.tmp_dir, "nxmate"), "w") as fw:
            for update_item in self.update_list:
                fw.write(update_item["infile1"] + "\n")
                fw.write(update_item["infile2"] + "\n")

    def get_command(self, cfg):
        command = [commands_parser.Command(STAGE=self.STAGE_NAME,
                                           path=sys.executable,
                                           args=[
                                               os.path.join(self.python_modules_home, "spades_pipeline", "scripts",
                                                            "preprocess_nxmate_reads.py"),
                                               "--args_filename", os.path.join(self.tmp_dir, "nxmate"),
                                               "--dst", self.dst,
                                               "--threads", str(options_storage.args.threads)],
                                           short_name=self.short_name)]
        return command

    def __init__(self, dir_for_split_reads, tmp_dir, *args):
        super(PreprocessNxmateReads, self).__init__(*args)
        self.tmp_dir = tmp_dir
        self.dst = dir_for_split_reads


class PreprocessContigs(stage.Stage):
    STAGE_NAME = "Preprocess additional contigs"

    # (gzipped, old_filename, new_filename)
    update_list = []

    def process_Ns_in_additional_contigs(self, dataset_data, dst, log):
        self.dst = dst
        for reads_library in dataset_data:
            if reads_library["type"].endswith("contigs"):
                new_entry = []
                for contigs in reads_library["single reads"]:
                    if contigs in options_storage.dict_of_prefixes:
                        ext = options_storage.dict_of_prefixes[contigs]
                        basename = contigs
                    else:
                        basename, ext = os.path.splitext(contigs)

                    gzipped = False
                    if ext.endswith(".gz"):
                        gzipped = True
                        if contigs not in options_storage.dict_of_prefixes:
                            basename, _ = os.path.splitext(basename)
                    new_filename = os.path.join(dst, os.path.basename(basename) + ".fasta")
                    if contigs in options_storage.dict_of_prefixes:
                        del options_storage.dict_of_prefixes[contigs]
                    new_entry.append(new_filename)
                    self.update_list.append({"gzipped": gzipped, "old_filename": contigs, "new_filename": new_filename})
                reads_library["single reads"] = new_entry

    def generate_config(self, cfg):
        self.process_Ns_in_additional_contigs(self.dataset_data, self.dst, self.log)
        with open(os.path.join(self.tmp_dir, "contigs"), "w") as fw:
            for update_item in self.update_list:
                fw.write(str(update_item["gzipped"]) + "\n")
                fw.write(update_item["old_filename"] + "\n")
                fw.write(update_item["new_filename"] + "\n")

    def get_command(self, cfg):
        command = [commands_parser.Command(STAGE=self.STAGE_NAME,
                                           path=sys.executable,
                                           args=[
                                               os.path.join(self.python_modules_home, "spades_pipeline", "scripts",
                                                            "preprocess_contigs.py"),
                                               "--args_filename", os.path.join(self.tmp_dir, "contigs"),
                                               "--dst", self.dst,
                                               "--threshold_for_breaking_additional_contigs",
                                               str(options_storage.THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS)],
                                           short_name=self.short_name)]
        return command

    def __init__(self, dir_for_split_reads, tmp_dir, *args):
        super(PreprocessContigs, self).__init__(*args)
        self.tmp_dir = tmp_dir
        self.dst = dir_for_split_reads


# splitting interlaced reads and processing Ns in additional contigs if needed
class PreprocessReadsStage(stage.Stage):
    STAGE_NAME = "Preprocess reads"
    stages = []

    def __init__(self, cfg, *args):
        super(PreprocessReadsStage, self).__init__(*args)

        self.dir_for_split_reads = os.path.join(options_storage.args.output_dir, "split_input")
        self.tmp_dir = os.path.join(self.dir_for_split_reads, "tmp")

        if support.dataset_has_interlaced_reads(self.dataset_data):
            self.stages.append(PreprocessInterlacedReads(self.dir_for_split_reads, self.tmp_dir, "preprocess_12",
                                                          self.output_files, self.tmp_configs_dir,
                                                          self.dataset_data, self.log,
                                                          self.bin_home,
                                                          self.ext_python_modules_home,
                                                          self.python_modules_home))

        if support.dataset_has_nxmate_reads(self.dataset_data):
            self.stages.append(PreprocessNxmateReads(self.dir_for_split_reads, self.tmp_dir, "preporocess_nxmate",
                                                     self.output_files, self.tmp_configs_dir,
                                                     self.dataset_data, self.log,
                                                     self.bin_home,
                                                     self.ext_python_modules_home,
                                                     self.python_modules_home))

        if support.dataset_has_additional_contigs(self.dataset_data):
            self.stages.append(PreprocessContigs(self.dir_for_split_reads, self.tmp_dir, "preprocess_ac",
                                                 self.output_files, self.tmp_configs_dir,
                                                 self.dataset_data, self.log,
                                                 self.bin_home,
                                                 self.ext_python_modules_home,
                                                 self.python_modules_home))

        options_storage.args.dataset_yaml_filename = os.path.join(options_storage.args.output_dir,
                                                                  "input_dataset.yaml")

        cfg["dataset"].yaml_filename = options_storage.args.dataset_yaml_filename

    def generate_config(self, cfg):
        if not os.path.isdir(self.dir_for_split_reads):
            os.makedirs(self.dir_for_split_reads)
        if not os.path.isdir(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        for stage in self.stages:
            stage.generate_config(cfg)

        if sys.version.startswith("2."):
            import pyyaml2 as pyyaml
        elif sys.version.startswith("3."):
            import pyyaml3 as pyyaml

        with open(options_storage.args.dataset_yaml_filename, 'w') as f:
            pyyaml.dump(self.dataset_data, f,
                        default_flow_style=False, default_style='"', width=float("inf"))

    def get_command(self, cfg):
        return [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name + "_start")] + \
               [x for stage in self.stages for x in stage.get_command(cfg)] + \
               [commands_parser.Command(STAGE=self.STAGE_NAME,
                                        path="true",
                                        args=[],
                                        short_name=self.short_name + "_finish",
                                        del_after=[os.path.relpath(self.tmp_dir, options_storage.args.output_dir)])]


def add_to_pipeline(pipeline, cfg, output_files, tmp_configs_dir, dataset_data, log,
                    bin_home, ext_python_modules_home, python_modules_home):
    if support.dataset_has_interlaced_reads(dataset_data) or support.dataset_has_additional_contigs(dataset_data) \
            or support.dataset_has_nxmate_reads(dataset_data):
        pipeline.add(PreprocessReadsStage(cfg, "preprocess", output_files, tmp_configs_dir,
                                          dataset_data, log, bin_home, ext_python_modules_home, python_modules_home))
