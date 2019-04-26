PathRacer: racing profile HMM paths on assembly graph
MANUAL

Overview
PathRacer is assembly graph against profile HMM aligning tool supporting
both local-local and global-local (aka glocal) alignment and both nucleotide and amino acid profile HMMs.
The tool finds all proper alignments rather than only the best one.
That allows extracting all genes satisfying HMM gene model from the assembly.



Input
For this moment the tool supports de Bruijn graphs in GFA format produced by SPAdes.
Contact us if you need some other format support.
Profile HMM should be in HMMER3 format, but one can pass nucleotide or ami acid sequence(s) to be converted to pHMM automatically.
k (de Bruijn overlap size) for the input graph also should be passed.

Output
For each pHMM (gene) the tool reports:
- <gene_name>.seqs.fa: sequences correspondent to K (parameter) best score paths ordered by score along with their alignment in CIGAR format
- <gene_name>.nucs.fa: (for amino acids pHHMs only) the same sequences in nucleotides
- <gene_name>.edges.fa: unique unitig (edge) paths correspondent to best score paths
- <gene_name>.{domtblout, pfamtblout, tblout}: (optional) unitig paths realignment by HMMER hmmalign in various formats
- event_graph_<gene_name>_component_<component_id>_size_<component_size>.cereal: (optional, debug output) connected components of the aligned graph
In addition:
- all.edges.fa: unique unitig paths for all pHMMs in one file
- pathracer.log: log file
- graph_with_hmm_paths.gfa: (optional) input graph with annotated unitig paths


Command line options



Examples




Restrictions and precautions
PathRacer is in intensive development. While we believe that the current release is reliable and stable,
some options could be changed or removed in further versions.

Running on large datasets requires a large amount of resources --- CPU, disk, memory, and stack.
Memory limit exceeding and stack overflow cause instant fall probably with no any informative message.
In case of that consider memory and stack limit increasing:
ulimit -s unlimited
export OMP_STACKSIZE=1G
./pathracer -m 500G ...
For extremely large and complicated datasets also consider running in a single thread:
./pathracer -t1 ... 

In case of any troubles, do not hesitate contact SPAdes support <> or the authors directly 
shlemovalex@gmail.com attaching the log file.
Your suggestions are also very welcome!

References:
If you are using PathRacer in your research, please refer to
bioARXIV










