#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os

# for restarting SPAdes
original_k_mers = None

dict_of_prefixes = dict()
dict_of_rel2abs = dict()

correct_scaffolds = False
run_truseq_postprocessing = False

SUPPORTED_PYTHON_VERSIONS = ["2.7", "3.2+"]  # major.minor format only, close ("-") and open ("+") ranges allowed
# allowed reads extensions for BayesHammer and for thw whole SPAdes pipeline
BH_ALLOWED_READS_EXTENSIONS = [".fq", ".fastq", ".bam", ".fq.gz", ".fastq.gz"]
IONTORRENT_ONLY_ALLOWED_READS_EXTENSIONS = [".bam"]
CONTIGS_ALLOWED_READS_EXTENSIONS = [".fa", ".fasta", ".fa.gz", ".fasta.gz"]
GRAPH_ALLOWED_READS_EXTENSIONS = [".gfa"]
ALLOWED_READS_EXTENSIONS = BH_ALLOWED_READS_EXTENSIONS + CONTIGS_ALLOWED_READS_EXTENSIONS + GRAPH_ALLOWED_READS_EXTENSIONS

# we support up to MAX_LIBS_NUMBER libs for each type of short-reads libs
MAX_LIBS_NUMBER = 9
OLD_STYLE_READS_OPTIONS = ["--12", "-1", "-2", "-s", "--merged"]
SHORT_READS_TYPES = {"pe": "paired-end", "s": "single", "mp": "mate-pairs", "hqmp": "hq-mate-pairs"}
# other libs types:
LONG_READS_TYPES = ["pacbio", "sanger", "nanopore", "trusted-contigs", "untrusted-contigs", "fl-rna"]
GRAPH_READS_TYPES = ["assembly-graph"]

SHORT_STAGES_NAME = ["ec", "as", "mc", "scc", "tpp"]

# final contigs and scaffolds names
contigs_name = "contigs.fasta"
scaffolds_name = "scaffolds.fasta"
secondary_scaffolds_name = "raw_scaffolds.fasta"
assembly_graph_name = "assembly_graph.fastg"
assembly_graph_name_gfa = "assembly_graph_with_scaffolds.gfa"
contigs_paths = "contigs.paths"
secondary_contigs_name = "raw_contigs.fasta"
scaffolds_paths = "scaffolds.paths"
secondary_scaffolds_paths = "raw_scaffolds.paths"
transcripts_name = "transcripts.fasta"
transcripts_paths = "transcripts.paths"
filtering_types = ["hard", "soft"]
bgc_stats_name = "hmm_statistics.txt"
gene_clusters_name = "gene_clusters.fasta"
domain_graph_name = "domain_graph.dot"

pipeline_state_dir = "pipeline_state"
biosyntheticspades_hmms = "biosynthetic_spades_hmms"
coronaspades_hmms = "coronaspades_hmms"

# other constants
MIN_K = 1
MAX_K = 127
RNA_MIN_K = 29
RNA_MAX_LOWER_K = 55
RNA_VIRAL_MAX_LOWER_K = 45
THRESHOLD_FOR_BREAKING_SCAFFOLDS = 3
THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS = 10
GAP_CLOSER_ENABLE_MIN_K = 55
SCC_K = 21

# default values constants
THREADS = 16
MEMORY = 250
K_MERS_RNA = [33, 49]
K_MERS_SHORT = [21, 33, 55]
K_MERS_150 = [21, 33, 55, 77]
K_MERS_250 = [21, 33, 55, 77, 99, 127]
K_MERS_PLASMID_100 = [21, 33, 55, 77]
K_MERS_PLASMID_LONG = [21, 33, 55, 77, 99, 127]

ITERATIONS = 1
TMP_DIR = "tmp"

READS_TYPES_USED_IN_CONSTRUCTION = ["paired-end", "single", "hq-mate-pairs"]
READS_TYPES_USED_IN_RNA_SEQ = ["paired-end", "single", "trusted-contigs", "untrusted-contigs", "pacbio", "nanopore", "fl-rna"]

BASE_STAGE = "read_conversion"
LAST_STAGE = "last"

first_command_line = None
args = None

original_dataset_data = None

# get path to checkpoint stage file
def get_stage_filename(stage_num, stage_short_name):
    stage_file_name = "stage_%d_%s" % (stage_num, stage_short_name)
    stage_checkpoint_path = os.path.join(args.output_dir, pipeline_state_dir, stage_file_name)
    return stage_checkpoint_path


# kmers were set by default, not SC, not IonTorrent data and not rna and temporary not meta (except metaplasmid)
def auto_K_allowed():
    return not args.k_mers and not args.single_cell and not args.iontorrent and not (args.meta and not args.plasmid)

def hmm_mode():
    return args.bio or args.custom_hmms or args.corona
