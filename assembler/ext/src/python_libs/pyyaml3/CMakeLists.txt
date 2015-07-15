############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(pyyaml3)

install(FILES __init__.py composer.py constructor.py
             cyaml.py dumper.py emitter.py  error.py
             events.py loader.py nodes.py parser.py  
             reader.py representer.py resolver.py  
             scanner.py serializer.py tokens.py
        DESTINATION share/spades/pyyaml3
        COMPONENT runtime)
