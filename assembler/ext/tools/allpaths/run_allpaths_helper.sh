#!/bin/sh

if [ $# -eq 0 ]
then
    echo "It it auxiliary script don't use it!"
    exit
fi

BUILD_DIR=${PWD}/../../../build/ext/allpaths/
WORK_DIR=${BUILD_DIR}/temp

APSRC=${BUILD_DIR}/allpaths_src    
WGSIMSRC=${BUILD_DIR}/wgsim_src 

if [ $1 = "make" ]
then    
    if [ ! -d ${BUILD_DIR} ]
    then 
        mkdir -p ${BUILD_DIR}
    fi
    if [ ! -d ${APSRC} ]
    then
        echo "copying allpaths"
        cp -r ./allpaths_src ${APSRC}
        ./${APSRC}/configure >${BUILD_DIR}/make.log 2>${BUILD_DIR}/make.err
        cp -r ./wgsim_src ${WGSIMSRC}
    fi

    echo "== making ALLPATHS =="
    cd ${APSRC}      
    make >>${BUILD_DIR}/make.log 2>>${BUILD_DIR}/make.err
    cd ${APSRC}/src/allpaths_cache/
    chmod +x *.pl   
    ln -s -t . $(find ../ -perm /u=x,g=x,o=x) 2> /dev/null
    cd ${WGSIMSRC}
    make >/dev/null
    echo "== making finished =="

    exit
fi

export PATH=$PATH:$APSRC/src/allpaths_cache/:$APSRC/src/
export PERL5LIB=$APSRC/src/allpaths_cache/

if [ $1 = "prepare" ]
then   
    echo "== preparing data =="
    sh ${WORK_DIR}/prepare.sh >${WORK_DIR}/prepare.out 2>${WORK_DIR}/prepare.err
    echo "== preparing data finished =="

    exit
fi

if [ $1 = "assemble" ]
then   
    echo "== assembling data =="
    sh ${WORK_DIR}/assemble.sh >${WORK_DIR}/assemble.out 2>${WORK_DIR}/assemble.err
    echo "== assembling data finished =="

    exit
fi
