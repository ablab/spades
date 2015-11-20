#!/bin/bash
virtualenv venv
. venv/bin/activate
pip install Flask
pip install Flask-Session
pip install jsonpickle
pip install pyparsing
