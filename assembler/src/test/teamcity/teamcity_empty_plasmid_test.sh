#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


if [[ -s $1 ]] ; then
echo "$1 was not empty!"
exit 43
else
echo "$1 is empty."
exit 0
fi ;
