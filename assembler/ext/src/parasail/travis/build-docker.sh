#!/bin/bash
set -e -x

SUFFIX=$1
echo "SUFFIX=$SUFFIX"

cd /io
name=parasail-`./version.sh`-manylinux1_$SUFFIX
autoreconf -fi
./configure --prefix=`pwd`/$name
make -j 2
make install
tar -cvzf $name.tar.gz $name
mkdir -p dist
cp $name.tar.gz dist
