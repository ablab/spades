############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades_pipeline)

# SPAdes pipeline scripts
install(FILES hammer_logic.py process_cfg.py spades_logic.py corrector_logic.py support.py options_storage.py lucigen_nxmate.py
        DESTINATION share/spades/spades_pipeline
        COMPONENT runtime)
# TruSpades module
install(FILES truspades/reference_construction.py truspades/moleculo_filter_contigs.py truspades/break_by_coverage.py truspades/moleculo_postprocessing.py truspades/barcode_extraction.py truspades/generate_quality.py truspades/id_generation.py truspades/launch_options.py truspades/string_dist_utils.py
        DESTINATION share/spades/spades_pipeline/truspades
        COMPONENT runtime)
# Common module
install(FILES common/alignment.py common/parallel_launcher.py common/sam_parser.py common/SeqIO.py common/__init__.py
        DESTINATION share/spades/spades_pipeline/common
        COMPONENT runtime)
