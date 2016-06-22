############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(hattrie C)

add_library(hattrie STATIC
            ahtable.c hat-trie.c misc.c murmurhash3.c)

