#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

current_folder=$1
etalons_folder=$2

echo "Regenerating etalon folder $etalons_folder from folder $current_folder"

if [ ! -d $current_folder ]; then
    echo "Error: source folder $current_folder not found"
    exit 9
fi

rm -rf $etalons_folder
mkdir $etalons_folder

for f in $current_folder/K*
    do
		target_folder=$etalons_folder/$(basename $f)
		echo "Copying saves from subfolder $(basename $f)"
        mkdir $target_folder
		cp -r $f/saves $target_folder
    done
echo "Done"
