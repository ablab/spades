#!/bin/bash

#***************************************************************************
##* Copyright (c) 2011-2014 Saint-Petersburg Academic University
##* All Rights Reserved
##* See file LICENSE for details.
#****************************************************************************


grep "SPAdes-$1.$2.$3-Darwin.tar.gz" spades$1$2$3.txt  > mac$1$2$3.txt
grep "SPAdes-$1.$2.$3-Linux.tar.gz" spades$1$2$3.txt  > linux$1$2$3.txt
grep "SPAdes-$1.$2.$3.tar.gz" spades$1$2$3.txt  > src$1$2$3.txt
grep "/node/362/done?sid=" spades$1$2$3.txt  > reg$1$2$3.txt
grep "/manual.html" spades$1$2$3.txt  > manual$1$2$3.txt
grep "curl" spades$1$2$3.txt  > curl$1$2$3.txt
grep "Wget" spades$1$2$3.txt  > wget$1$2$3.txt
echo "Total users" 
./filter.py spades$1$2$3.txt
echo "Sources" 
./filter.py src$1$2$3.txt
echo "Linux" 
./filter.py linux$1$2$3.txt
echo "Mac OS" 
./filter.py mac$1$2$3.txt
echo "Registration form" 
./filter.py reg$1$2$3.txt
echo "Link from manual" 
./filter.py manual$1$2$3.txt
echo "curl" 
./filter.py curl$1$2$3.txt
echo "wget" 
./filter.py wget$1$2$3.txt

