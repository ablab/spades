# HMM-guided mode
The majority of SPAdes assembly modes (normal multicell, single-cell, rnaviral, meta and of course biosynthetic) also supports HMM-guided mode as implemented in biosyntheticSPAdes. The detailed description could be found in [biosyntheticSPAdes paper](https://genome.cshlp.org/content/early/2019/06/03/gr.243477.118), but in short: amino acid profile HMMs are aligned to the edges of assembly graph. After this the subgraphs containing the set of matches ("domains") are extracted and all possible paths through the domains that are supported both by paired-end data (via scaffolds) and graph topology are obtained (putative biosynthetic gene clusters).

HMM-guided mode could be enabled via providing a set of HMMs via `--custom-hmms` option. In HMM guided mode the set of contigs and scaffolds (see [SPAdes output](output.md#spades-output) section for more information ) is kept intact, however additional [biosyntheticSPAdes output](output.md#biosyntheticspades-output) represents the output of HMM-guided assembly.

Note that normal biosyntheticSPAdes mode (via `--bio` option) is a bit different from HMM-guided mode: besides using the special set of profile HMMS representing a family of NRSP/PKS domains also includes a set of assembly graph simplification and processing settings aimed for fuller recovery of biosynthetic gene clusters.

## coronaSPAdes mode

Given an increased interest in coronavirus research we developed a coronavirus assembly mode for SPAdes assembler (also known as coronaSPAdes). It allows to assemble full-length coronaviridae genomes from the transcriptomic and metatranscriptomic data. Algorithmically, coronaSPAdes is an rnaviralSPAdes that uses the set of HMMs from [Pfam SARS-CoV-2 2.0](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam_SARS-CoV-2_2.0/) set as well as additional HMMs as outlined by [(Phan et al, 2019)](https://doi.org/10.1093/ve/vey035). coronaSPAdes could be run via a dedicated `coronaspades.py` script. See [coronaSPAdes preprint](https://www.biorxiv.org/content/10.1101/2020.07.28.224584v1) for more information about rnaviralSPAdes,  coronaSPAdes and HMM-guided mode. Output for any HMM-related mode (`--bio`, `--corona`, or `--custom-hmms` flags) is the same with biosyntheticSPAdes' output.


## wastewaterSPAdes mode

SARS-CoV-2 wastewater samples are extensively collected and studied because it allows quantitative assessment of viral load in surrounding populations. We developed wastewaterSPAdes that solves SARS-CoV-2 deconvolution problem using assembly graph structure.
To use wastewaterSPAdes, you'll need to:

- Set `--sewage` flag to the `coronaspades.py`.
- Provide the SARS-CoV-2 reference genome as trusted contigs.

Results of wastewaterSPAdes are stored in `lineages.csv` file. First column contains strain name, and second column contains estimated abundance of this strain in the sample.
