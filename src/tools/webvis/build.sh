
############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#!/bin/bash
virtualenv venv
. venv/bin/activate
pip install Flask
pip install Flask-Session
pip install jsonpickle
pip install pyparsing
