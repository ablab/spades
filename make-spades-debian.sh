rm -r spades-*
rm spades_*

VERSION="$(cat assembler/VERSION)"
mkdir spades-$VERSION
cp -r assembler/src spades-$VERSION/src
cp -r assembler/configs spades-$VERSION/configs
cp -r assembler/debian spades-$VERSION/debian
cp -r assembler/ext spades-$VERSION/ext
cp -r assembler/test_dataset spades-$VERSION/test_dataset
cp assembler/LICENSE spades-$VERSION
cp assembler/README spades-$VERSION
cp assembler/VERSION spades-$VERSION
cp assembler/makefileDebian spades-$VERSION/makefile
cp assembler/spades.py spades-$VERSION
cp assembler/spades_config.info.template spades-$VERSION
cp assembler/spades_download_binary.py spades-$VERSION
cp assembler/spades_download_bayeshammer.py spades-$VERSION
cp assembler/spades_init.py spades-$VERSION

cd spades-$VERSION
rm src/CMakeListsInternal.txt
touch src/CMakeListsInternal.txt
rm -r src/abruijn
rm -r src/test
rm -r configs/debruijn/datasets_archive
rm -r configs/debruijn/datasets
find . -name ".?*" | xargs rm -r
rm -r src/tools/blast-fasta
rm -r src/tools/clean_contaminations
rm -r src/tools/consensus
rm -r src/tools/contig_analysis
rm -r src/tools/cuckoo_test
rm -r src/tools/divide_fastq
rm -r src/tools/estimation
rm -r src/tools/estimation_deprecated
rm -r src/tools/filter
rm -r src/tools/kmerstats
rm -r src/tools/mapreads
rm -r src/tools/maps_test
rm -r src/tools/quake-correct
rm -r src/tools/ukonnen
rm -r src/tools/copy_ungzip.sh
rm -r src/tools/dot_conv.sh
rm -r src/tools/run_velvet_unpaired.sh
rm -r src/tools/run_velvet2.sh
rm -r src/tools/quality/libs/mauve
rm -r src/tools/quality/libs/gage
rm -r src/tools/quality/libs/genemark_suite_linux_64
rm -r src/tools/quality/libs/report
rm -r src/tools/spades_build/

cp configs/debruijn/config.info.template configs/debruijn/config.info
cp configs/hammer/config.info.template configs/hammer/config.info
cp configs/debruijn/detail_info_printer.info.template configs/debruijn/detail_info_printer.info
cp configs/debruijn/simplification.info.template configs/debruijn/simplification.info
cp configs/debruijn/distance_estimation.info.template configs/debruijn/distance_estimation.info
cp configs/debruijn/long_contigs/lc_config.info.template configs/debruijn/long_contigs/lc_config.info
cp configs/debruijn/long_contigs/lc_params.info.template configs/debruijn/long_contigs/lc_params.info
cp spades_config.info.template spades_config.info

debuild -us -uc

cd ..

