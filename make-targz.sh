rm -r spades-*
rm spades_*

VERSION="$(cat assembler/VERSION)"
mkdir -p spades-$VERSION/src/tools

cp -r assembler/src/debruijn spades-$VERSION/src/debruijn
cp -r assembler/src/hammer spades-$VERSION/src/hammer
cp -r assembler/src/include spades-$VERSION/src/include
cp -r assembler/src/io spades-$VERSION/src/io
cp -r assembler/src/tools/quality spades-$VERSION/src/tools/quality
cp -r assembler/src/tools/spades_pipeline spades-$VERSION/src/tools/spades_pipeline
cp assembler/src/CMakeLists.txt spades-$VERSION/src/CMakeLists.txt

cp -r assembler/configs spades-$VERSION/configs
cp -r assembler/debian spades-$VERSION/debian
cp -r assembler/ext spades-$VERSION/ext
cp -r assembler/test_dataset spades-$VERSION/test_dataset
cp assembler/LICENSE spades-$VERSION
cp assembler/README spades-$VERSION
cp assembler/VERSION spades-$VERSION
cp assembler/makefileDebian spades-$VERSION/makefile
cp assembler/spades.py spades-$VERSION
cp assembler/quast.py spades-$VERSION
cp assembler/spades_config.info.template spades-$VERSION
cp assembler/spades_download_binary.py spades-$VERSION
cp assembler/spades_download_bayeshammer.py spades-$VERSION
cp assembler/spades_init.py spades-$VERSION
cp assembler/manual.html spades-$VERSION
cp assembler/quality.html spades-$VERSION

cd spades-$VERSION
touch src/CMakeListsInternal.txt
rm -r configs/debruijn/datasets_archive
rm -r configs/debruijn/datasets
find . -name ".?*" | xargs rm -r

cp configs/debruijn/config.info.template configs/debruijn/config.info
cp configs/hammer/config.info.template configs/hammer/config.info
cp configs/debruijn/detail_info_printer.info.template configs/debruijn/detail_info_printer.info
cp configs/debruijn/simplification.info.template configs/debruijn/simplification.info
cp configs/debruijn/distance_estimation.info.template configs/debruijn/distance_estimation.info
cp configs/debruijn/long_contigs/lc_config.info.template configs/debruijn/long_contigs/lc_config.info
cp configs/debruijn/long_contigs/lc_params.info.template configs/debruijn/long_contigs/lc_params.info
cp spades_config.info.template spades_config.info
mv src/include/k_range.hpp.template src/include/k_range.hpp

cd ..

tar -pczf spades-$VERSION.tar.gz spades-$VERSION
