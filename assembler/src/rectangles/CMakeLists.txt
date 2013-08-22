############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(rectangles)

install(FILES N50.py bigraph.py diagonal.py fastaparser.py not_used_code.py
              rectangle.py test_util.py abstract_graph.py check_diags.py  
              experimental.py graph.py pathsets.py rectangle_set.py  
              saveparser.py utils.py
        DESTINATION share/spades/rectangles
        COMPONENT runtime)
install(FILES rrr.py
        DESTINATION share/spades/rectangles
        COMPONENT runtime)
