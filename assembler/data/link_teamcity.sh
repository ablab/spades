#!/bin/bash
set -e
if [ ! -e input ];
then
    ln -s ../../../../input input
fi