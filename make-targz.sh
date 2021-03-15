VERSION="$(cat assembler/VERSION)"
TARGET_DIR=SPAdes-$VERSION
rm -rf $TARGET_DIR
SRC_DIR=$TARGET_DIR/src
mkdir -p $SRC_DIR

cp -r assembler/src/common $SRC_DIR/
cp -r assembler/src/projects $SRC_DIR/
cp -r assembler/src/include $SRC_DIR/
cp -r assembler/src/cmake $SRC_DIR/
cp -r assembler/src/spades_pipeline $SRC_DIR/
cp assembler/src/CMakeLists.txt $SRC_DIR/

cp -r assembler/configs $TARGET_DIR/configs
cp -r assembler/ext $TARGET_DIR/ext

# cleaning .pyc and .pyo
rm -f $SRC_DIR/*.pyc
rm -f $SRC_DIR/*.pyo
rm -rf $SRC_DIR/__pycache__/
rm -f $SRC_DIR/spades_pipeline/*.pyc
rm -f $SRC_DIR/spades_pipeline/*.pyo
rm -rf $SRC_DIR/spades_pipeline/__pycache__/
rm -f $SRC_DIR/spades_pipeline/*/*.pyo
rm -f $SRC_DIR/spades_pipeline/*/*.pyc
rm -rf $SRC_DIR/spades_pipeline/*/__pycache__/
rm -f $TARGET_DIR/ext/src/python_libs/*/*.pyc
rm -f $TARGET_DIR/ext/src/python_libs/*/*.pyo
rm -rf $TARGET_DIR/ext/src/python_libs/*/__pycache__/

cp -r assembler/biosynthetic_spades_hmms $TARGET_DIR/
cp -r assembler/coronaspades_hmms $TARGET_DIR/
cp -r assembler/test_dataset $TARGET_DIR/test_dataset
cp -r assembler/test_dataset_plasmid $TARGET_DIR/test_dataset_plasmid
cp assembler/LICENSE $TARGET_DIR/
cp README.md $TARGET_DIR/
cp assembler/VERSION $TARGET_DIR/
cp assembler/spades.py $TARGET_DIR/
cp assembler/rnaviralspades.py $TARGET_DIR/
cp assembler/rnaspades.py $TARGET_DIR/
cp assembler/metaspades.py $TARGET_DIR/
cp assembler/plasmidspades.py $TARGET_DIR/
cp assembler/metaviralspades.py $TARGET_DIR/
cp assembler/metaplasmidspades.py $TARGET_DIR/
cp assembler/coronaspades.py $TARGET_DIR/
cp assembler/spades_compile.sh $TARGET_DIR/
cp assembler/spades_init.py $TARGET_DIR/
cp assembler/manual.html $TARGET_DIR/
cp assembler/rnaspades_manual.html $TARGET_DIR/
cp assembler/changelog.html $TARGET_DIR/
cp assembler/GPLv2.txt $TARGET_DIR/

cd $TARGET_DIR/
touch src/CMakeListsInternal.txt
find . -name ".?*" | xargs rm -r

cd ..

tar -pczf SPAdes-$VERSION.tar.gz SPAdes-$VERSION
rm -r SPAdes-$VERSION
