#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os


class OptionStorage:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if hasattr(self, '_initialized'):
            # Prevent re-initialization
            return
        self._initialized = True
        # for restarting SPAdes
        self.original_k_mers = None

        self.dict_of_prefixes = dict()
        self.dict_of_rel2abs = dict()

        self.MINIMAL_PYTHON_VERSION = (3, 8)
        # allowed reads extensions for BayesHammer and for thw whole SPAdes pipeline
        self.BH_ALLOWED_READS_EXTENSIONS = [".fq", ".fastq", ".bam", ".fq.gz", ".fastq.gz"]
        self.IONTORRENT_ONLY_ALLOWED_READS_EXTENSIONS = [".bam"]
        self.CONTIGS_ALLOWED_READS_EXTENSIONS = [".fa", ".fasta", ".fa.gz", ".fasta.gz"]
        self.GRAPH_ALLOWED_READS_EXTENSIONS = [".gfa"]
        self.ALLOWED_READS_EXTENSIONS = (self.BH_ALLOWED_READS_EXTENSIONS + self.CONTIGS_ALLOWED_READS_EXTENSIONS +
                                         self.GRAPH_ALLOWED_READS_EXTENSIONS + [".sra"])

        # we support up to MAX_LIBS_NUMBER libs for each type of short-reads libs
        self.MAX_LIBS_NUMBER = 9
        self.OLD_STYLE_READS_OPTIONS = ["--12", "-1", "-2", "-s", "--merged"]
        self.SHORT_READS_TYPES = {"pe": "paired-end", "s": "single", "mp": "mate-pairs", "hqmp": "hq-mate-pairs"}
        # other libs types:
        self.LONG_READS_TYPES = ["pacbio", "sanger", "nanopore", "trusted-contigs", "untrusted-contigs", "fl-rna"]
        self.GRAPH_READS_TYPES = ["assembly-graph"]

        self.SHORT_STAGES_NAME = ["ec", "as", "mc", "scc", "tpp"]

        # final contigs and scaffolds names
        self.contigs_name = "contigs.fasta"
        self.scaffolds_name = "scaffolds.fasta"
        self.secondary_scaffolds_name = "raw_scaffolds.fasta"
        self.assembly_graph_name = "assembly_graph.fastg"
        self.assembly_graph_name_gfa = "assembly_graph_with_scaffolds.gfa"
        self.contigs_paths = "contigs.paths"
        self.secondary_contigs_name = "raw_contigs.fasta"
        self.scaffolds_paths = "scaffolds.paths"
        self.secondary_scaffolds_paths = "raw_scaffolds.paths"
        self.transcripts_name = "transcripts.fasta"
        self.transcripts_paths = "transcripts.paths"
        self.filtering_types = ["hard", "soft"]
        self.bgc_stats_name = "hmm_statistics.txt"
        self.gene_clusters_name = "gene_clusters.fasta"
        self.domain_graph_name = "domain_graph.dot"
        self.sewage_lineages = "lineages.csv"

        self.pipeline_state_dir = "pipeline_state"
        self.biosyntheticspades_hmms = "biosynthetic_spades_hmms"
        self.coronaspades_hmms = "coronaspades_hmms"

        # other constants
        self.MIN_K = 1
        self.MAX_K = 127
        self.RNA_MIN_K = 29
        self.RNA_MAX_LOWER_K = 55
        self.RNA_VIRAL_MAX_LOWER_K = 45
        self.THRESHOLD_FOR_BREAKING_SCAFFOLDS = 3
        self.THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS = 10
        self.GAP_CLOSER_ENABLE_MIN_K = 55
        self.SCC_K = 21

        # default values constants
        self.THREADS = 16
        self.MEMORY = 250
        self.K_MERS_RNA = [33, 49]
        self.K_MERS_SHORT = [21, 33, 55]
        self.K_MERS_150 = [21, 33, 55, 77]
        self.K_MERS_250 = [21, 33, 55, 77, 99, 127]
        self.K_MERS_PLASMID_100 = [21, 33, 55, 77]
        self.K_MERS_PLASMID_LONG = [21, 33, 55, 77, 99, 127]

        self.ITERATIONS = 1
        self.TMP_DIR = "tmp"

        self.READS_TYPES_USED_IN_CONSTRUCTION = ["paired-end", "single", "hq-mate-pairs"]
        self.READS_TYPES_USED_IN_RNA_SEQ = ["paired-end", "single", "trusted-contigs", "untrusted-contigs", "pacbio", "nanopore", "fl-rna"]

        self.BASE_STAGE = "read_conversion"
        self.LAST_STAGE = "last"

        self.first_command_line = None
        self.args = None
        self.spades_version = ""

        self.original_dataset_data = None

    # get path to checkpoint stage file
    def get_stage_filename(self, stage_num, stage_short_name):
        stage_file_name = "stage_%d_%s" % (stage_num, stage_short_name)
        stage_checkpoint_path = os.path.join(self.args.output_dir, self.pipeline_state_dir, stage_file_name)
        if not os.path.exists(os.path.dirname(stage_checkpoint_path)):
            os.makedirs(os.path.dirname(stage_checkpoint_path))
        return stage_checkpoint_path

    # kmers were set by default, not SC, not IonTorrent data and not rna and temporary not meta (except metaplasmid)
    def auto_K_allowed(self):
        return (not self.args.k_mers and not self.args.single_cell and
                not self.args.iontorrent and not (self.args.meta and not self.args.plasmid))

    def hmm_mode(self):
        return self.args.bio or self.args.custom_hmms or self.args.corona
