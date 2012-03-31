OLDEXT=dot
NEWEXT=${1/#.}

if [ $# -ge 2 ]
then
	cd $2	
fi

find . -iname "*.${OLDEXT}" |
while read F
do
	NEWFILE="${F/%${OLDEXT}/${NEWEXT}}"
	echo converting $NEWFILE ...
	dot "${F}" -T${NEWEXT} > "${NEWFILE}"
done
