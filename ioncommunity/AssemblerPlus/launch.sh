#!/bin/bash
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved
#AUTORUNDISABLED

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
 
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
 
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

VERSION="5.0.0" # major.minor.bug

# set the ram
#$ -l mem_free=${PLUGINCONFIG__RAM},s_vmem=${PLUGINCONFIG__RAM}

# ===================================================
# Plugin functions
# ===================================================

#*! @function
#  @param  $*  the command to be executed
run ()
{
	echo "running: $*";
	eval $*;
	EXIT_CODE="$?";
}

#*! @function
set_output_paths ()
{
	PLUGIN_OUT_BAM_NAME=`echo ${TSP_FILEPATH_BAM} | sed -e 's_.*/__g' | sed -e 's/bam/basecaller.bam/g'`; 
}

# ===================================================
# Plugin initialization
# ===================================================

# Set defaults
set_output_paths;
ASSEMBLER_PATH="${DIRNAME}/bin/";

# remove some files if they are there
run "rm -rf ${TSP_FILEPATH_PLUGIN_DIR}/*.html ${TSP_FILEPATH_PLUGIN_DIR}/info*.json"

# ===================================================
# Run AssemblerPlus Plugin
# ===================================================
#barcoded run
if [ -f ${TSP_FILEPATH_BARCODE_TXT} ]; then
    echo "${TSP_FILEPATH_BARCODE_TXT}";

    BARCODE_LINES=`wc -l ${TSP_FILEPATH_BARCODE_TXT} | awk '{print \$1}'`;
    echo ${BARCODE_LINES};
    NUM_BARCODES=`expr ${BARCODE_LINES} - 2`;

#    echo "${NUM_BARCODES} barcodes found. Files with fewer than ${PLUGINCONFIG__MIN_READS} reads will be skipped." > ${TSP_FILEPATH_PLUGIN_DIR}/samples_block.html;
    
    CTR=0;
    IFS=$'\n';
    for BARCODE_LINE in `cat ${TSP_FILEPATH_BARCODE_TXT} | grep "^barcode"`
    do
        BARCODE_ID=`echo ${BARCODE_LINE} | awk 'BEGIN{FS=","} {print $2}'`;
        if [ -n "$PLUGINCONFIG__ONLY_BARCODES" ] && [[ ,$PLUGINCONFIG__ONLY_BARCODES, != *,$BARCODE_ID,* ]]; then
          continue
        fi

        BARCODE_SEQ=`echo ${BARCODE_LINE} | awk 'BEGIN{FS=","} {print $3}'`;
        BARCODE_BAM_NAME="${BARCODE_ID}_${PLUGIN_OUT_BAM_NAME}";

	echo "";
	echo "";
	echo "#########################";
	echo "#Starting assembly of ${BARCODE_BAM_NAME}";
	echo "#########################";
	echo "";
	echo "";
	
	#see if bam file exists
	if [ -f ${BASECALLER_DIR}/${BARCODE_BAM_NAME} ]; then
	    #create sub dir
	    if [ ! -f ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ} ]; then
		run "mkdir -p ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}";
	    fi

	    #remove any existing bam links that may be there
	    if [ -f ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME} ]; then
		run "rm ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME}";
	    fi

	    #create sym link to bam file
	    if [ ! -f ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME} ]; then
		run "ln -snf ${BASECALLER_DIR}/${BARCODE_BAM_NAME} ${TSP_FILEPATH_PLUGIN_DIR}/${BARCODE_ID}.${BARCODE_SEQ}/${BARCODE_BAM_NAME}";
	    fi

            #build call to the assembler.pl script which will take care of the rest
	    run "python ${ASSEMBLER_PATH}/RunAssembler.py \"${BARCODE_ID}\" \"${BARCODE_SEQ}\" \"${BARCODE_BAM_NAME}\"";
	    CTR=`expr ${CTR} + 1`;
	fi
    done

    run "python ${ASSEMBLER_PATH}/GenerateReport.py ${TSP_FILEPATH_PLUGIN_DIR}/info*.json"

#nonbarcoded run
else
    echo "";
    echo "";
    echo "#########################";
    echo "#Starting assembly of ${PLUGIN_OUT_BAM_NAME}";
    echo "#########################";
    echo "";
    echo "";

    #remove any link that may be there
    if [ -f ${TSP_FILEPATH_PLUGIN_DIR}/${PLUGIN_OUT_BAM_NAME} ]; then
	run "rm ${TSP_FILEPATH_PLUGIN_DIR}/${PLUGIN_OUT_BAM_NAME}";
    fi

    #create the sym link
    if [ ! -f ${TSP_FILEPATH_PLUGIN_DIR}/${PLUGIN_OUT_BAM_NAME} ]; then
	run "ln -snf ${BASECALLER_DIR}/${PLUGIN_OUT_BAM_NAME} ${TSP_FILEPATH_PLUGIN_DIR}/${PLUGIN_OUT_BAM_NAME}";
    fi

    run "python ${ASSEMBLER_PATH}/RunAssembler.py ${PLUGIN_OUT_BAM_NAME}"
    run "python ${ASSEMBLER_PATH}/GenerateReport.py ${TSP_FILEPATH_PLUGIN_DIR}/info.json"
fi
