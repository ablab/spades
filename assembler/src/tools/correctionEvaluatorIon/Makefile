all: bam2fastq comparator cluster_summaries error_sites fasta_slice singletons ion_dist_table

DC=dmd
INCLUDES=-IBioD
D_OPTS=-O -release -inline
RDMD_OPTS=--compiler=$(DC) --build-only

.PHONY: clean

clean:
	rm -f bam2fastq comparator cluster_summaries error_sites fasta_slice singletons ion_dist_table

bam2fastq: bam2fastq.d
	rdmd $(RDMD_OPTS) $(D_OPTS) $(INCLUDES) bam2fastq.d

comparator: comparator.d
	rdmd $(RDMD_OPTS) $(D_OPTS) $(INCLUDES) comparator.d

cluster_summaries: cluster_summaries.d
	rdmd $(RDMD_OPTS) $(D_OPTS) $(INCLUDES) cluster_summaries.d

error_sites: error_sites.d
	rdmd $(RDMD_OPTS) $(D_OPTS) $(INCLUDES) error_sites.d

fasta_slice: fasta_slice.d
	rdmd $(RDMD_OPTS) $(D_OPTS) $(INCLUDES) fasta_slice.d

singletons: singletons.d
	rdmd $(RDMD_OPTS) $(D_OPTS) $(INCLUDES) singletons.d

ion_dist_table: ion_dist_table.d
	rdmd $(RDMD_OPTS) $(D_OPTS) $(INCLUDES) ion_dist_table.d
