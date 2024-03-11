#!/bin/bash

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -p|--port)
    PORT="$2"
    shift
    ;;
    -c|--clean)
    CLEAN=YES
    ;;
    *)
    # unknown option
    ;;
esac
shift
done

if [ $CLEAN ];
then
  echo Cleaning leftovers...
  rm -rf session/*
  rm -f flask_session/*
fi

. venv/bin/activate
python index.py $PORT
