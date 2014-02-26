#!/bin/bash

#***************************************************************************
##* Copyright (c) 2011-2014 Saint-Petersburg Academic University
##* All Rights Reserved
##* See file LICENSE for details.
#****************************************************************************


grep "SPAdes-[0-9].[0-9].[0-9]-Darwin.tar.gz" $1 > mac.txt
grep "SPAdes-[0-9].[0-9].[0-9]-Linux.tar.gz" $1  > linux.txt
grep "SPAdes-[0-9].[0-9].[0-9].tar.gz" $1  > src.txt
grep "/node/362/done?sid=" $1  > reg.txt
grep "/manual.html" $1  > manual.txt
grep "curl" $1  > curl.txt
grep "Wget" $1  > wget.txt
echo "Total users" 
./filter.py $1
echo "Sources" 
./filter.py src.txt
echo "Linux" 
./filter.py linux.txt
echo "Mac OS" 
./filter.py mac.txt
echo "Registration form" 
./filter.py reg.txt
echo "Link from manual" 
./filter.py manual.txt
echo "curl" 
./filter.py curl.txt
echo "wget" 
./filter.py wget.txt

