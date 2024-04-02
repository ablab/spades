# Installation


SPAdes requires a 64-bit Linux system or Mac OS and Python (supported versions are Python 2.7, and Python3: 3.2 and higher) to be pre-installed on it. To obtain SPAdes you can either download binaries or download source code and compile it yourself. []()

In case of successful installation the following files will be placed in the `bin` directory:

-   `spades.py` (main executable script)
-   `metaspades.py` (main executable script for [metaSPAdes](#meta))
-   `plasmidspades.py` (main executable script for [plasmidSPAdes](#plasmid))
-   `metaplasmidspades.py` (main executable script for [metaplasmidSPAdes](#metaextrachromosomal))
-   `metaviralspades.py` (main executable script for [metaviralSPAdes](#metaextrachromosomal))
-   `rnaspades.py` (main executable script for [rnaSPAdes](rnaspades_manual.html))
-   `rnaviralspades.py` (main executable script for rnaviralSPAdes)
-   `coronaspades.py` (wrapper script for coronaSPAdes mode)
-   `spades-core`  (assembly module)
-   `spades-gbuilder`  (standalone graph builder application)
-   `spades-gmapper`  (standalone long read to graph aligner)
-   `spades-kmercount`  (standalone k-mer counting application)
-   `spades-hammer`  (read error correcting module for Illumina reads)
-   `spades-ionhammer`  (read error correcting module for IonTorrent reads)
-   `spades-bwa`  ([BWA](http://bio-bwa.sourceforge.net) alignment module which is required for mismatch correction)
-   `spades-corrector-core`  (mismatch correction module)

<a name="sec2.1"></a>
## Downloading SPAdes Linux binaries

To download [SPAdes Linux binaries](http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz) and extract them, go to the directory in which you wish SPAdes to be installed and run:

``` bash

    wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz
    tar -xzf SPAdes-3.15.4-Linux.tar.gz
    cd SPAdes-3.15.4-Linux/bin/
```

In this case you do not need to run any installation scripts - SPAdes is ready to use. We also suggest adding SPAdes installation directory to the `PATH` variable. []()

Note, that pre-build binaries do not work on new Linux kernels.

<a name="sec2.2"></a>
## Downloading SPAdes binaries for Mac

To obtain [SPAdes binaries for Mac](http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Darwin.tar.gz), go to the directory in which you wish SPAdes to be installed and run:

``` bash

    curl http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Darwin.tar.gz -o SPAdes-3.15.4-Darwin.tar.gz
    tar -zxf SPAdes-3.15.4-Darwin.tar.gz
    cd SPAdes-3.15.4-Darwin/bin/
```

Just as in Linux, SPAdes is ready to use and no further installation steps are required. We also suggest adding SPAdes installation directory to the `PATH` variable. []()

<a name="sec2.3"></a>
## Downloading and compiling SPAdes source code

If you wish to compile SPAdes by yourself you will need the following libraries to be pre-installed:

-   g++ (version 5.3.1 or higher)
-   cmake (version 3.5 or higher)
-   zlib
-   libbz2

If you meet these requirements, you can download the [SPAdes source code](http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4.tar.gz):

``` bash

    wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4.tar.gz
    tar -xzf SPAdes-3.15.4.tar.gz
    cd SPAdes-3.15.4
```

and build it with the following script:

``` bash

    ./spades_compile.sh
```

SPAdes will be built in the directory `./bin`. If you wish to install SPAdes into another directory, you can specify full path of destination folder by running the following command in `bash` or `sh`:

``` bash

    PREFIX=<destination_dir> ./spades_compile.sh
```

for example:

``` bash

    PREFIX=/usr/local ./spades_compile.sh
```

which will install SPAdes into `/usr/local/bin`.

After installation you will get the same files (listed above) in `./bin` directory (or `<destination_dir>/bin` if you specified PREFIX). We also suggest adding `bin` directory to the `PATH` variable. []()

<a name="sec2.4"></a>
## Verifying your installation

For testing purposes, SPAdes comes with a toy data set (reads that align to first 1000 bp of *E. coli*). To try SPAdes on this data set, run:

``` bash

    <spades installation dir>/bin/spades.py --test
```

If you added `bin` folder from SPAdes installation directory to the `PATH` variable, you can run:

``` bash

    spades.py --test
```

For the simplicity we further assume that `bin` folder from SPAdes installation directory is added to the `PATH` variable.

If the installation is successful, you will find the following information at the end of the log:

``` plain

===== Assembling finished. Used k-mer sizes: 21, 33, 55

 * Corrected reads are in spades_test/corrected/
 * Assembled contigs are in spades_test/contigs.fasta
 * Assembled scaffolds are in spades_test/scaffolds.fasta
 * Assembly graph is in spades_test/assembly_graph.fastg
 * Assembly graph in GFA format is in spades_test/assembly_graph_with_scaffolds.gfa
 * Paths in the assembly graph corresponding to the contigs are in spades_test/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in spades_test/scaffolds.paths

======= SPAdes pipeline finished.

========= TEST PASSED CORRECTLY.

SPAdes log can be found here: spades_test/spades.log

Thank you for using SPAdes!
```
