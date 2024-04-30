# HMM-guided mode
The majority of SPAdes assembly modes (multicell, single-cell, rnaviral, meta and biosynthetic) also supports HMM-guided mode as implemented in biosyntheticSPAdes. 
The detailed description could be found in [biosyntheticSPAdes paper](https://genome.cshlp.org/content/early/2019/06/03/gr.243477.118), but in short: amino acid profile HMMs are aligned to the edges of assembly graph. 
After this the subgraphs containing the set of matches ("domains") are extracted and all possible paths through the domains that are supported both by paired-end data (via scaffolds) and graph topology are obtained (putative biosynthetic gene clusters).

HMM-guided mode is enabled via providing a set of HMMs (`*.hmm.gz` file) via `--custom-hmms` option. In HMM-guided mode the set of contigs and scaffolds (see [SPAdes output](output.md#spades-output) section for more information ) is kept intact, however additional [biosyntheticSPAdes output](output.md#biosyntheticspades-output) represents the output of HMM-guided assembly.

We provide an example of HMM utility in viral assembly, along with general advice on constructing HMM profile sets for various purposes, in our [paper on noroviral assembly](https://www.mdpi.com/2079-7737/12/8/1066).

Note that normal biosyntheticSPAdes mode (via `--bio` option) is a bit different from HMM-guided mode: besides using the special set of profile HMMS representing a family of NRSP/PKS domains also includes a set of assembly graph simplification and processing settings aimed for fuller recovery of biosynthetic gene clusters.


## coronaSPAdes mode

Given an increased interest in coronavirus research we developed a coronavirus assembly mode for SPAdes assembler (also known as coronaSPAdes). 
It allows to assemble full-length coronaviridae genomes from the transcriptomic and metatranscriptomic data. Algorithmically, coronaSPAdes is an rnaviralSPAdes that uses the set of HMMs from [Pfam SARS-CoV-2 2.0](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam_SARS-CoV-2_2.0/) set as well as additional HMMs as outlined by [(Phan et al, 2019)](https://doi.org/10.1093/ve/vey035). coronaSPAdes could be run via a dedicated `coronaspades.py` script. See [coronaSPAdes paper](https://academic.oup.com/bioinformatics/article/38/1/1/6354349) for more information about rnaviralSPAdes, coronaSPAdes and HMM-guided mode. Output for any HMM-related mode (`--bio`, `--corona`, or `--custom-hmms` flags) is the same with biosyntheticSPAdes' output.



