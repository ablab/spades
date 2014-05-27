############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades_pipeline)

# SPAdes pipeline scripts
install(FILES hammer_logic.py process_cfg.py spades_logic.py dipspades_logic.py support.py options_storage.py
        DESTINATION share/spades/spades_pipeline
        COMPONENT runtime)
# Mismatch corrector module
install(FILES corrector.py
        DESTINATION share/spades/spades_pipeline
        COMPONENT runtime)

