############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(joblib3)

install(FILES __init__.py disk.py format_stack.py
             func_inspect.py _compat.py hashing.py
             logger.py  memory.py my_exceptions.py 
             numpy_pickle.py  parallel.py testing.py
        DESTINATION share/spades/joblib3
        COMPONENT runtime)
