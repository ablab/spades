############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades_pipeline)

# SPAdes pipeline scripts
install(FILES
        run_contig_breaker.py
        commands_parser.py
        process_cfg.py
        support.py
        options_storage.py
        options_parser.py
        DESTINATION share/spades/spades_pipeline
        COMPONENT runtime)

install(DIRECTORY scripts
        DESTINATION share/spades/spades_pipeline/
        COMPONENT runtime
        PATTERN "*.pyc" EXCLUDE)


install(DIRECTORY stages
        DESTINATION share/spades/spades_pipeline/
        COMPONENT runtime
        PATTERN "*.pyc" EXCLUDE)

# TruSpades module
install(DIRECTORY truspades
        DESTINATION share/spades/spades_pipeline/
        COMPONENT runtime
        PATTERN "*.pyc" EXCLUDE)

# Common module
install(DIRECTORY common
        DESTINATION share/spades/spades_pipeline/
        COMPONENT runtime
        PATTERN "*.pyc" EXCLUDE)

# Executors scripts
install(DIRECTORY executors
        DESTINATION share/spades/spades_pipeline/
        COMPONENT runtime
        PATTERN "*.pyc" EXCLUDE)
