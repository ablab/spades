############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

OLDEXT=dot
NEWEXT=${1/#.}

if [ $# -ge 2 ]
then
	if [ -d "$2" ]
	then
		cd $2	
	else
		echo "Couldn't find folder $2. Exiting."
		exit
	fi
fi

find . -iname "*.${OLDEXT}" |
while read F
do
	NEWFILE="${F/%${OLDEXT}/${NEWEXT}}"
	echo converting $NEWFILE ...
	dot "${F}" -T${NEWEXT} > "${NEWFILE}"
done
