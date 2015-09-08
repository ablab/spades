#!/bin/bash

if [ -n $1 ];
then
  pkill online_vis
  rm flask_session/*
fi

. venv/bin/activate
python index.py
