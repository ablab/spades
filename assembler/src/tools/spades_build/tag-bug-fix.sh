#!/bin/bash
set -e

#eval "git pull --rebase"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FILE=$DIR/bug-fix.version

if [ -z "$1" ];
then 
    read bug_fix_version < $FILE
else 
    $bug_fix_version=$1
fi

bug_fix_version=$(($bug_fix_version+1))
echo "Bug fixed number is set to " $bug_fix_version
echo $bug_fix_version > $FILE



eval "git stash"
eval "git add "$FILE
eval "git ci -m'"


