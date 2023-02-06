############################################################################
# Copyright (c) 2023 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

install(PROGRAMS "spades.py"
        DESTINATION bin
        COMPONENT spades)
install(PROGRAMS "metaspades.py"
        DESTINATION bin
        COMPONENT spades)
install(PROGRAMS "plasmidspades.py"
        DESTINATION bin
        COMPONENT spades)
install(PROGRAMS "rnaspades.py"
        DESTINATION bin
        COMPONENT spades)
install(PROGRAMS "metaviralspades.py"
        DESTINATION bin
        COMPONENT spades)
install(PROGRAMS "metaplasmidspades.py"
        DESTINATION bin
        COMPONENT spades)
install(PROGRAMS "rnaviralspades.py"
        DESTINATION bin
        COMPONENT spades)
install(PROGRAMS "coronaspades.py"
        DESTINATION bin
        COMPONENT spades)

install(FILES "spades_init.py"
        DESTINATION bin
        COMPONENT spades)

# SPAdes pipeline scripts
install(DIRECTORY "spades_pipeline"
        DESTINATION share/spades
        COMPONENT runtime
        PATTERN "*.pyc" EXCLUDE
        PATTERN "__pycache__" EXCLUDE)
