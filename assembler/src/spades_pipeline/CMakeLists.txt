############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades_pipeline)

# SPAdes pipeline scripts
install(FILES bh_aux.py bh_logic.py process_cfg.py spades_logic.py support.py options_storage.py
        DESTINATION share/spades/spades_pipeline
        COMPONENT runtime)
# Mismatch corrector module
install(FILES corrector.py
        DESTINATION share/spades/spades_pipeline
        COMPONENT runtime
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ)

