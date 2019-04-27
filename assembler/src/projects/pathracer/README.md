PathRacer: racing profile HMM paths on assembly graph
=====================================================
MANUAL
------

### Overview
<!-- PathRacer is assembly graph against profile HMM aligning tool supporting -->
<!-- both _local-local_ and _global-local_ (aka _glocal_) alignment and both nucleotide and amino acid profile HMMs. -->
<!-- The tool finds all proper alignments rather than only the best one. -->
<!-- That allows extracting all genes satisfying HMM gene model from the assembly. -->
<!--  -->
**PathRacer** is a novel standalone tool that aligns profile HMM directly to the
assembly graph (performing the codon translation on fly for amino acid pHMMs).
The tool provides the set of most probable paths traversed by a HMM through the
whole assembly graph, regardless whether the sequence of interested is encoded
on the single contig or scattered across the set of edges, therefore
significantly improving the recovery of sequences of interest even from
fragmented metagenome assemblies.

### Input
For this moment the tool supports only _de Bruijn_ graphs in GFA format produced by **SPAdes**.
Contact us if you need some other format support.

Profile HMM should be in **HMMer3** format, but one can pass nucleotide or amino acid sequence(s) to be converted to pHMM(s) automatically.

### Output
For each pHMM (gene) the tool reports:

- **&lt;gene\_name&gt;.seqs.fa**: sequences correspondent to _N_ (parameter, see below) best score paths ordered by score along with their alignment in CIGAR format
- **&lt;gene\_name&gt;.nucs.fa**: _(for amino acids pHHMs only)_ the same sequences in nucleotides
- **&lt;gene\_name&gt;.edges.fa**: unique unitig (edge) paths correspondent to best score paths above
- **&lt;gene\_name&gt;.{domtblout, pfamtblout, tblout}**: _(optional)_ unitig paths realignment by **HMMer3** `hmmalign` in various formats
- **event\_graph\_&lt;gene\_name&gt;\_component\_&lt;component\_id&gt;\_size\_&lt;component\_size&gt;.cereal**: _(optional, debug output)_ connected components of the aligned graph
- **&lt;component\_id&gt;.dot**: _(optional, plot)_ connected component of matched neighborhood subgraph
- **&lt;component\_id&gt;\_&lt;path\_index&gt;.dot**: _(optional, plot)_ neighborhood of the found path

In addition:

- **all.edges.fa**: unique unitig paths for all pHMMs in one file
- **pathracer.log**: log file
- **graph\_with\_hmm\_paths.gfa**: _(optional)_ input graph with annotated unitig paths


### Command line options
Required positional arguments:

1. Query gene models file (.hmm file or .fasta)
2. Graph in GFA format
3. _k_ (_de Bruijn_ overlap size) for the input graph

Main options:

- `--output`, `-o` DIR: output directory
- `--hmm` | `--nt` | `--aa`: match against pHMM(s) [default] | nucleotide sequences | amino acid sequences
- `--queries` Q1 [Q2 [...]]: queries names to lookup [default: all queries from input query file]
- `--global` | `--local`: perform HMM-global, graph-local (aka _glocal_, default) or HMM-local, graph-local HMM matching
- `--length`, `-l` L: minimal length of resultant matched sequence; if &le;1 then to be multiplied on aligned HMM length [default: 0.9]
- `--top` N: extract up to _N_ top paths [default: 10000]
- `--rescore`: rescore paths by **HMMer3**
- `--threads`, `-t` T: the number of parallel threads [default: 16]
- `--parallel-components`: process connected components of neighborhood subgraph in parallel
- `--memory`, `-m` M: RAM limit for PathRacer in GB (terminates if exceeded) [default: 100]

Debug output control:

- `--debug`: enable extensive debug console output
- `--draw`: draw pictures around the interesting edges
- `--annotate-graph`: emit paths in GFA graph
- `--export-event-graph`: export Event Graph in .cereal format

Heuristics options:

- `--max-size` MAX\_SIZE: maximal component size to consider [default: INF]
- `--max-insertion-length`: maximal allowed number of successive I-emissions [default: 30]
- `--no-top-score-filter`: disable top score Event Graph vertices filter. Increases sensitivity of deep analysis (`--top` &gt; 50000)

_In addition:_ Some other developer options that are not supposed to be tuned by end-user. Could be removed in further releases.

### Examples
One can download example datasets from here <http://cab.spbu.ru/software/pathracer/>

- **urban_strain.gfa**: strain assembly graph of Singapore clinical isolation ward wastewater metagenome (SRA accession SRR5997548, dataset H1)
- **bla\_all.hmm**: 159 beta-lactamase family pHMMs from **AMRFinder** database <https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/>
- **ecoli\_mc.gfa**: assembly graph of _E. coli str. K12 substr. MG1655_ multicell isolate (SRA accession ERA000206)
- **bac.hmm**: _16S_, _23S_, and _5S_ ribosomal subunit pHMMs
- **synth\_strain\_gbuilder.gfa**: strain assembly graph of SYNTH mock metagenome dataset (SRA accession SRR606249)
- **synth16S\_new.fa**: _16S_ rRNA sequences of SYNTH organisms obtained from [SILVA](https://www.arb-silva.de/) database

Lookup for beta-lactamase genes (amino acid pHMMs) in Singapore wastewater  
```
./pathracer bla_all.hmm urban_strain.gfa 55 --output pathracer_urban_strain_bla_all
```

Lookup for _16S_/_5S_/_23S_ (nucleotide HMMs) in _E.coli_ multicell assembly  
```
./pathracer bac.hmm ecoli_mc.gfa 55 --output pathracer_ecoli_mc_bac
```

Look up for known _16S_ sequences in _E.coli_ multicell assembly  
```
./pathracer synth16S_new.fa ecoli_mc.gfa 55 --nt --output pathracer_ecoli_mc_16S_seqs
```

Look up for known _16S_ sequences in SYNTH mock metagenome assembly  
```
./pathracer synth16S_new.fa synth_strain_gbuilder.gfa 55 --nt --output pathracer_synth_strain_gbuider_16S_seqs
```

Let us extract **all** _16S_ sequences from SYNTH mock metagenome assembly.
For this we increase `--top` and disable Event Graph vertices filter (`--no-top-score-filter`)
Deep analysis of extremely complicated dataset also require stack and memory limits tuning  
```
ulimit -s unlimited  
export OMP_STACKSIZE=1G  
./pathracer bac.hmm synth_strain_gbuilder.gfa 55 --queries 16S_rRNA -m 250 --top 1000000 --output pathracer_synth_strain_gbuilder_16s --no-top-score-filter
```

### References
If you are using **PathRacer** in your research, please refer to <https://www.biorxiv.org/content/10.1101/562579v1>

In case of any problems running **PathRacer** please contact SPAdes support <spades.support@cab.spbu.ru> attaching the log file.
Your suggestions are also very welcome!











