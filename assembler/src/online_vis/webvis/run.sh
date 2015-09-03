#!/bin/bash

if [ $1 == clean ];
then
  pkill online_vis
  rm flask_session/*
fi

. venv/bin/activate
python index.py
