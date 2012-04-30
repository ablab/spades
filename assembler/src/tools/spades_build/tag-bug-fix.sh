#!/bin/bash
set -e -x

VERSION="$(cat assembler/VERSION)"
BRANCH=${VERSION%.*}

eval "git checkout spades_$BRANCH"

echo "BRANCH set to " $BRANCH

eval "git pull --rebase"
eval "git stash"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FILE=$DIR/bug-fix.version

if [ -z "$1" ];
then
    read FIX < $FILE
else
    $FIX=$1
fi

FIX=$(($FIX+1))

cd assembler
dch -v $BRANCH.$FIX
cd ..

echo $FIX > $FILE
echo $BRANCH.$FIX > assembler/VERSION


git diff assembler/VERSION $FILE assembler/debian/changelog | less

git commit assembler/VERSION $FILE assembler/debian/changelog

cat assembler/VERSION

eval "git push"

eval "git tag spades_$BRANCH.$FIX"

eval "git push --tag"

eval "git stash pop"