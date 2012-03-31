OLDEXT=dot
NEWEXT=${1/#.}

find . -iname "*.${OLDEXT}" |
while read F
do
	NEWFILE="${F/%${OLDEXT}/${NEWEXT}}"
	echo converting $NEWFILE ...
	dot "${F}" -T${NEWEXT} > "${NEWFILE}"
done
