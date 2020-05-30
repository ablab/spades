This is a support software release for metaviralSPAdes publication

It contains script for viral assembly from metagenomes (assembler/metaviralspades.py), which is based on metaplasmidSPAdes.

For metaviral release no binary distribution is available now - you should install from source code. 

-   g++ (version 5.3.1 or higher)
-   cmake (version 2.8.12 or higher)
-   zlib
-   libbz2

should be installed before metaviralSPAdes.

To install you'll need to run

``` bash

    ./assembler/spades_compile.sh
```

metaviralSPAdes binaries will be built in the directory `./bin`.


Input/output options  are same as in regular SPAdes -- see http://cab.spbu.ru/files/release3.14.1/manual.html for details.


Additional output (results of repeat resolution without any bacterial contigs removed) is available as <Output_folder>/<MAX_K>/before_chromosome_removal.fasta

