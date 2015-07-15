#!/bin/bash

set -e
if [ ! -e input ];
then
    ln -s /tmp/data/input input
fi
