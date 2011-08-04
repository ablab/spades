#!/bin/bash

#renaming templates to originals
cp -n .project.template .project
cp -n .cproject.template .cproject

cd src/debruijn
cp -n .project.template .project
cp -n .cproject.template .cproject

cd ../abruijn
cp -n .project.template .project
cp -n .cproject.template .cproject

cd ../hammer
cp -n .project.template .project
cp -n .cproject.template .cproject

# going to root folder (assembler)
cd ../..
./prepare_debug.sh
./prepare_release.sh
