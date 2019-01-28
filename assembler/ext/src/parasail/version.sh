#!/bin/sh
#
# This script extracts the version from parasail.h, which is the master
# location for this information.
#
if [ ! -f parasail.h ]; then
    echo "version.sh: error: parasail.h does not exist" 1>&2
    exit 1
fi
MAJOR=`egrep '^#define +PARASAIL_VERSION_MAJOR +[0-9]+$' parasail.h`
MINOR=`egrep '^#define +PARASAIL_VERSION_MINOR +[0-9]+$' parasail.h`
PATCH=`egrep '^#define +PARASAIL_VERSION_PATCH +[0-9]+$' parasail.h`
if [ -z "$MAJOR" -o -z "$MINOR" -o -z "$PATCH" ]; then
    echo "version.sh: error: could not extract version from parasail.h" 1>&2
    exit 1
fi
MAJOR=`echo $MAJOR | awk '{ print $3 }'`
MINOR=`echo $MINOR | awk '{ print $3 }'`
PATCH=`echo $PATCH | awk '{ print $3 }'`
echo $MAJOR.$MINOR.$PATCH | tr -d '\n'

