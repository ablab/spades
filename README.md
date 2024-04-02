# About SPAdes

SPAdes is an assembly toolkit containing various assembly pipelines.

- [Complete SPAdes user manual](https://ablab.github.io/spades/)

- [SPAdes download page](https://github.com/ablab/spades/releases/)

- [Latest SPAdes publication](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102)


# Quick start

- Complete user manual can be found [here](https://ablab.github.io/spades/). Information below is provided merely for your convenience and cannot be considered as the user guide.

- SPAdes is an assembler for second-generation sequencing data (Illumina or IonTorrent). PacBio and Nanopore reads are supported *only* as supplementary data. SPAdes can assemble genomes, metagenomes, transcriptomes, viral genomes etc. 

- Download SPAdes binaries for [Linux](https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz) or [MacOS](https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Darwin.tar.gz). You can also compile SPAdes from [source](https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5.tar.gz) (requires g++ 9.0+, cmake 3.16+, zlib and libbz2). SPAdes requires only Python 3.8+ to be installed.

- Test your SPAdes installation by running

```
    bin/spades.py --test
```

- A single paired-end library (separate files, gzipped):

```
    bin/spades.py -1 left.fastq.gz -2 right.fastq.gz -o output_folder
```

- A single paired-end library (interlaced reads):

```
    bin/spades.py --12 interlaced.fastq -o output_folder
```

- Two paired-end libraries (separate files):

```
    bin/spades.py --pe1-1 1_left.fastq --pe1-2 1_right.fastq --pe2-1 2_left.fastq --pe2-2 2_right.fastq -o output_folder
```

- IonTorrent data:
```
    bin/spades.py --iontorrent -s it_reads.fastq -o output_folder
```

- A paired-end library coupled with long PacBio reads:

```
    bin/spades.py -1 left.fastq.gz -2 right.fastq.gz --pacbio pb.fastq -o output_folder
```

- Available assembly modes: `--isolate`, `--sc`, `--plasmid`, `--meta`, `--metaplasmid`, `--metaviral`, `--rna`, `--rnaviral`, `--corona`, `--bio`.


# Citation
If you use SPAdes in your research, please cite [our latest paper](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102).

In case you perform hybrid assembly using  PacBio or Nanopore reads, you may also cite [Antipov et al., 2015](http://bioinformatics.oxfordjournals.org/content/early/2015/11/20/bioinformatics.btv688.short). If you use multiple paired-end and/or mate-pair libraries you may additionally cite papers describing SPAdes repeat resolution algorithms [Prjibelski et al., 2014](http://bioinformatics.oxfordjournals.org/content/30/12/i293.short) and [Vasilinetc et al., 2015](http://bioinformatics.oxfordjournals.org/content/31/20/3262.abstract). 

If you use other pipelines, please cite the following papers:

-   metaSPAdes: [Nurk et al., 2017](https://genome.cshlp.org/content/27/5/824.short).
-   plasmidSPAdes: [Antipov et al., 2016](https://academic.oup.com/bioinformatics/article/32/22/3380/2525610).
-   metaplasmidSPAdes / plasmidVerify: [Antipov et al., 2019](https://genome.cshlp.org/content/29/6/961.short)
-   metaviralSPAdes / viralVerify: [Antipov et al., 2020](https://academic.oup.com/bioinformatics/article-abstract/36/14/4126/5837667)
-   rnaSPAdes: [Bushmanova et al., 2019](https://academic.oup.com/gigascience/article/8/9/giz100/5559527).
-   biosyntheticSPAdes: [Meleshko et al., 2019](https://genome.cshlp.org/content/early/2019/06/03/gr.243477.118?top=1).
-   coronaSPAdes paper is currently available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.28.224584v1.abstract).

You may also include older papers [Nurk, Bankevich et al., 2013](http://link.springer.com/chapter/10.1007%2F978-3-642-37195-0_13) or [Bankevich, Nurk et al., 2012](http://online.liebertpub.com/doi/abs/10.1089/cmb.2012.0021), especially if you assemble single-cell data.


# Feedback and bug reports

Please, leave your comments and bug reports at [our GitHub repository tracker](https://github.com/ablab/spades/issues). If you have any troubles running SPAdes, please attach us `params.txt` and `spades.log` from the output folder.

