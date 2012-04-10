rm -r spades-*
rm spades_*
cp -r assembler spades-2.0.0

cd spades-2.0.0
rm -r src/abruijn
rm -r src/bayes
rm -r src/test
rm -r data
rm -r configs/debruijn/datasets_archive
rm -r configs/debruijn/datasets
rm -r docs
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

rm gen_k
rm editcfg
rm install_prerequirements.sh
rm prepare_cfg
rm prepare_cfg_eclipse
rm run

mv makefileDebian makefile

debuild -us -uc

cd ..
scp spades_2.0.0* 192.168.222.223:spades2/

