VERSION="$(cat assembler/VERSION)"
rm -rf SPAdes-$VERSION
mkdir -p SPAdes-$VERSION/src

cp -r assembler/src/debruijn SPAdes-$VERSION/src/debruijn
cp -r assembler/src/hammer SPAdes-$VERSION/src/hammer
cp -r assembler/src/dipspades SPAdes-$VERSION/src/dipspades
cp -r assembler/src/ionhammer SPAdes-$VERSION/src/ionhammer
cp -r assembler/src/include SPAdes-$VERSION/src/include
cp -r assembler/src/io SPAdes-$VERSION/src/io
cp -r assembler/src/mph_index SPAdes-$VERSION/src/mph_index
cp -r assembler/src/cmake SPAdes-$VERSION/src/cmake
cp -r assembler/src/ssw SPAdes-$VERSION/src/ssw
cp -r assembler/src/spades_pipeline SPAdes-$VERSION/src/spades_pipeline
cp assembler/src/CMakeLists.txt SPAdes-$VERSION/src/CMakeLists.txt

cp -r assembler/configs SPAdes-$VERSION/configs
cp -r assembler/ext SPAdes-$VERSION/ext
rm -r SPAdes-$VERSION/ext/include/cute
rm -r SPAdes-$VERSION/ext/include/teamcity_boost

# cleaning .pyc and .pyo
rm -f SPAdes-$VERSION/src/spades_pipeline/*.pyc
rm -f SPAdes-$VERSION/src/spades_pipeline/*.pyo
rm -f SPAdes-$VERSION/ext/include/python_libs/joblib/*.pyc
rm -f SPAdes-$VERSION/ext/include/python_libs/joblib/*.pyo

cp -r assembler/test_dataset SPAdes-$VERSION/test_dataset
cp assembler/LICENSE SPAdes-$VERSION
cp assembler/README SPAdes-$VERSION
cp assembler/VERSION SPAdes-$VERSION
cp assembler/spades.py SPAdes-$VERSION
cp assembler/spades_compile.sh SPAdes-$VERSION
cp assembler/spades_init.py SPAdes-$VERSION
cp assembler/manual.html SPAdes-$VERSION
cp assembler/changelog.html SPAdes-$VERSION
cp assembler/GPLv2.txt SPAdes-$VERSION

cd SPAdes-$VERSION
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
cp configs/debruijn/log.properties.template configs/debruijn/log.properties
cp configs/debruijn/path_extend/pe_params.info.template configs/debruijn/path_extend/pe_params.info
cp configs/debruijn/coverage_based_rr.info.template   configs/debruijn/coverage_based_rr.info
cp configs/dipspades/config.info.template configs/dipspades/config.i
cp configs/dipspades/log.properties.template configs/dipspades/log.propertie
cp configs/ionhammer/ionhammer.cfg.template configs/ionhammer/ionhammer.cfg

cd ..

tar -pczf SPAdes-$VERSION.tar.gz SPAdes-$VERSION
rm -r SPAdes-$VERSION
