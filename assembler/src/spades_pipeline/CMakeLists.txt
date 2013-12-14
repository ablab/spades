############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades_pipeline)

# SPAdes pipeline scripts
install(FILES hammer_logic.py process_cfg.py spades_logic.py dipspades_logic.py corrector_logic.py support.py options_storage.py lucigen_nxmate.py
        DESTINATION share/spades/spades_pipeline
        COMPONENT runtime)
# Moleculo module
install(FILES alignment.py moleculo_filter_contigs.py break_by_coverage.py moleculo_postprocessing.py
        DESTINATION share/spades/spades_pipeline
        COMPONENT runtime)

