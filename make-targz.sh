
############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

VERSION="$(cat assembler/VERSION)"
: "${TARGET_DIR:=SPAdes-$VERSION}"
rm -rf $TARGET_DIR
SRC_DIR=$TARGET_DIR/src
mkdir -p $SRC_DIR

cp -r src/common $SRC_DIR/
cp -r src/projects $SRC_DIR/
cp -r src/include $SRC_DIR/
cp -r src/cmake $SRC_DIR/
cp -r src/spades_pipeline $SRC_DIR/
cp src/CMakeLists.txt $SRC_DIR/

cp -r configs $TARGET_DIR/configs
cp -r ext $TARGET_DIR/ext

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

cp -r biosynthetic_spades_hmms $TARGET_DIR/
cp -r coronaspades_hmms $TARGET_DIR/
cp -r test_dataset $TARGET_DIR/test_dataset
cp -r test_dataset_plasmid $TARGET_DIR/test_dataset_plasmid
cp LICENSE $TARGET_DIR/
cp README.md $TARGET_DIR/
cp VERSION $TARGET_DIR/
cp spades.py $TARGET_DIR/
cp rnaviralspades.py $TARGET_DIR/
cp rnaspades.py $TARGET_DIR/
cp metaspades.py $TARGET_DIR/
cp plasmidspades.py $TARGET_DIR/
cp metaviralspades.py $TARGET_DIR/
cp metaplasmidspades.py $TARGET_DIR/
cp coronaspades.py $TARGET_DIR/
cp spades_compile.sh $TARGET_DIR/
cp spades_init.py $TARGET_DIR/
cp changelog.html $TARGET_DIR/
cp GPLv2.txt $TARGET_DIR/

cd $TARGET_DIR
touch src/CMakeListsInternal.txt
find . -name ".?*" | xargs rm -r

cd ..

tar -pczf $TARGET_DIR.tar.gz $TARGET_DIR
rm -r $TARGET_DIR
