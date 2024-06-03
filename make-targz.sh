
############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

VERSION="$(cat VERSION)"
: "${TARGET_DIR:=SPAdes-$VERSION}"
rm -rf $TARGET_DIR
SRC_DIR=$TARGET_DIR/src
mkdir -p $SRC_DIR

cp -r src/common $SRC_DIR/
cp -r src/projects $SRC_DIR/
cp -r src/include $SRC_DIR/
cp -r src/cmake $SRC_DIR/
cp -r src/test $SRC_DIR/
cp src/CMakeLists.txt $SRC_DIR/

cp -r ext $TARGET_DIR/ext
cp -r docs $TARGET_DIR/docs

# cleaning .pyc and .pyo
shopt -s globstar
rm -f $SRC_DIR/**/*.pyc
rm -f $SRC_DIR/**/*.pyo
rm -fr $SRC_DIR/**/__pycache__/

cp LICENSE $TARGET_DIR/
cp README.md $TARGET_DIR/
cp VERSION $TARGET_DIR/
cp spades_compile.sh $TARGET_DIR/
cp changelog.md $TARGET_DIR/
cp GPLv2.txt $TARGET_DIR/

cd $TARGET_DIR
touch src/CMakeListsInternal.txt
find . -name ".?*" | xargs rm -r

cd ..

tar -pczf $TARGET_DIR.tar.gz $TARGET_DIR
rm -r $TARGET_DIR
