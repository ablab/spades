#!/bin/bash

#GLOBAL=$1

if hash easy_install 2>/dev/null; then
    sudo easy_install virtualenv
else
    sudo pip install virtualenv
fi
virtualenv venv
. venv/bin/activate
pip install Flask
pip install Flask-Session
pip install jsonpickle
pip install pyparsing
