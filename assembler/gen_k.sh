#!/bin/bash
FILENAME=src/debruijn/k.hpp

if [ $# = 0 ]
then 
echo "ERROR: K value should be specified as first parameter" 
else 
echo -e "#pragma once\n\nnamespace debruijn_graph {\n\tconst size_t K = "$1";\n}" > $FILENAME
fi
