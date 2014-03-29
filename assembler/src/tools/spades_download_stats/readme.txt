//***************************************************************************
///* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


To obtain statistics for SPAdes version X.Y.Z (starting 2.3.0).:

1. Run 
$> sudo ./grep_stat.sh X Y Z 
at morality.
For older versions run corresponding scripts.
You will obtain spadesXYZ.txt

2. Copy spadesXYZ.txt into directory with stat.sh and run
$> ./stat.sh X Y Z
or
$> ./stat_all.sh <path to your spadesXYZ.txt>
Unique users will be given in brackets.

If you wish to obtain statistics for several versions (e.g. the total number of unique users is not the sum of numbers of unique users for different versions) concatenate several files like spadesXYZ.txt into a single file and run:
$> ./stat_all.sh <path to the concatenated file>
