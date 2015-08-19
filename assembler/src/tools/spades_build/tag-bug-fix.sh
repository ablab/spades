#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e

VERSION="$(cat assembler/VERSION)"
BRANCH=${VERSION%.*}

eval "git checkout spades_$BRANCH"

echo "BRANCH set to " $BRANCH

eval "git pull --rebase"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FILE=$DIR/bug-fix.version


if [ -z "$1" ];
then
    FIX=`cat $FILE`
else
	$FIX=$1
fi

echo  $FIX 
cd assembler
dch -v $BRANCH.$FIX
cd ..

echo $BRANCH.$FIX > assembler/VERSION
FIX_NEXT=$(($FIX+1))
echo $FIX_NEXT > $FILE

git diff assembler/VERSION $FILE assembler/debian/changelog | less

git commit assembler/VERSION $FILE assembler/debian/changelog

cat assembler/VERSION

eval "git push"

eval "git tag spades_$BRANCH.$FIX"

eval "git push --tag"