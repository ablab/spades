rm -r spades-*

VERSION="$(cat assembler/VERSION)"
mkdir -p spades-$VERSION/src/tools

cp -r assembler/src/debruijn spades-$VERSION/src/debruijn
cp -r assembler/src/hammer spades-$VERSION/src/hammer
cp -r assembler/src/include spades-$VERSION/src/include
cp -r assembler/src/io spades-$VERSION/src/io
cp -r assembler/src/mph_index spades-$VERSION/src/mph_index
cp -r assembler/src/cmake spades-$VERSION/src/cmake
cp -r assembler/src/spades_pipeline spades-$VERSION/src/spades_pipeline
cp assembler/src/CMakeLists.txt spades-$VERSION/src/CMakeLists.txt

cp -r assembler/configs spades-$VERSION/configs
cp -r assembler/ext spades-$VERSION/ext
rm spades-$VERSION/ext/prepare_ext.sh
rm -r spades-$VERSION/ext/tools
rm -r spades-$VERSION/ext/src
rm -r spades-$VERSION/ext/include/cute
rm -r spades-$VERSION/ext/include/teamcity_boost

cp -r assembler/test_dataset spades-$VERSION/test_dataset
cp assembler/LICENSE spades-$VERSION
cp assembler/README spades-$VERSION
cp assembler/VERSION spades-$VERSION
cp assembler/spades.py spades-$VERSION
cp assembler/spades_config.info.template spades-$VERSION
cp assembler/spades_download_binary.py spades-$VERSION
cp assembler/spades_compile.sh spades-$VERSION
cp assembler/spades_init.py spades-$VERSION
cp assembler/manual.html spades-$VERSION

cd spades-$VERSION
touch src/CMakeListsInternal.txt
rm -r configs/debruijn/datasets_archive
rm -r configs/debruijn/datasets
rm  configs/debruijn/datasets*
rm  configs/debruijn/deprecated*
find . -name ".?*" | xargs rm -r

cp configs/debruijn/config.info.template configs/debruijn/config.info
cp configs/hammer/config.info.template configs/hammer/config.info
cp configs/debruijn/detail_info_printer.info.template configs/debruijn/detail_info_printer.info
cp configs/debruijn/simplification.info.template configs/debruijn/simplification.info
cp configs/debruijn/distance_estimation.info.template configs/debruijn/distance_estimation.info
cp configs/debruijn/path_extend/pe_params.info.template configs/debruijn/path_extend/pe_params.info

cd ..

tar -pczf spades-$VERSION.tar.gz spades-$VERSION
