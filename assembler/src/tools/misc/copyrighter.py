#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

##### for updating copyrights (e.g. from 2014 to 2015):  #####
### grep -rl "# Copyright (c) 2011-2014 Saint-Petersburg Academic University" . | xargs sed -i 's/# Copyright (c) 2011-2014 Saint-Petersburg Academic University/# Copyright (c) 2011-2015 Saint-Petersburg Academic University/g'
#####

##### for removing copyright (script-style, only SPbAU) #####
###  perl -0777 -i -pe 's/############################################################################\n# Copyright \(c\) 2011-2014 Saint-Petersburg Academic University\n# All Rights Reserved\n# See file LICENSE for details.\n############################################################################\n\n//igs' *.{py,sh}
#####

##### for removing copyright (code-style, two universities) #####
###  grep -rl "State University" . | grep -v "kmer_index.hpp" | xargs perl -0777 -i -pe 's/\Q\/\/***************************************************************************\E\n\Q\/\/* Copyright (c) 2015 Saint Petersburg State University\E\n\Q\/\/* Copyright (c) 2011-2014 Saint Petersburg Academic University\E\n\Q\/\/* All Rights Reserved\E\n\Q\/\/* See file LICENSE for details.\E\n\Q\/\/***************************************************************************\E\n\n//igs'
#####

##### for updating one-university copyright to two-university one #####
### grep -rl "Copyright (c) 2011-2014" . | xargs perl -0777 -i -pe 's/\Q\/\/***************************************************************************\E\n\Q\/\/* Copyright (c) 2011-2014 Saint-Petersburg Academic University\E\n\Q\/\/* All Rights Reserved\E\n\Q\/\/* See file LICENSE for details.\E\n\Q\/\/****************************************************************************\E\n/\/\/***************************************************************************\n\/\/* Copyright \(c\) 2015 Saint Petersburg State University\n\/\/* Copyright \(c\) 2011-2014 Saint Petersburg Academic University\n\/\/* All Rights Reserved\n\/\/* See file LICENSE for details.\n\/\/***************************************************************************\n/igs'
#####

##### for removing SPbAU copyrights from new files (specified in <LIST OF FILES>)
### for i in `cat <LIST OF FILES>`; do sed -i '' '/Saint Petersburg Academic University/d' $i; done
#####


import os
import shutil
import sys

script_comment = [
    '############################################################################',
    '# Copyright (c) 2015 Saint Petersburg State University',
    #'# Copyright (c) 2011-2014 Saint Petersburg Academic University',  # new copyrights should be with SPbSU only
    '# All Rights Reserved',
    '# See file LICENSE for details.',
    '############################################################################', 
    '']

code_comment = [
    '//***************************************************************************',
    '//* Copyright (c) 2015 Saint Petersburg State University',
    #'//* Copyright (c) 2011-2014 Saint Petersburg Academic University',  # new copyrights should be with SPbSU only
    '//* All Rights Reserved',
    '//* See file LICENSE for details.',
    '//***************************************************************************',
    '']

only_show = False


def insert_in_script(filename):
    print(filename)
    if only_show:
       return

    lines = open(filename).readlines()
    if (script_comment[1] + '\n') in lines:
        return

    modified = open(filename, 'w')
    if len(lines) and lines[0].startswith('#!'):
        modified.write(lines[0])
        modified.write('\n')
        lines = lines[1:]
    
    for com_line in script_comment:
        modified.write(com_line + '\n')
    
    for line in lines:
        modified.write(line)

    modified.close() 


def insert_in_code(filename):       
    print(filename)
    if only_show:
       return

    lines = open(filename).readlines()
    if (code_comment[1] + '\n') in lines:
        return

    modified = open(filename, 'w')
        
    for com_line in code_comment:
        modified.write(com_line + '\n')
    
    for line in lines:
        modified.write(line)

    modified.close()


def visit(arg, dirname, names):
    for name in names:
        path = os.path.join(dirname, name)
        if not os.path.isfile(path):
            continue
        if os.system('grep -Gq "Copyright (c) \([0-9]\{4\}-\)\{0,1\}[0-9]\{4\} Saint Petersburg" ' + path) == 0: # already copyrighted
            continue
        ext = os.path.splitext(name)[1]
        if arg and ext != arg:
            continue
        if (ext in ['.py', '.sh']) or name.lower().startswith('cmake'):
            insert_in_script(path)
        elif ext in ['.hpp', '.cpp', '.h', '.c']:
            insert_in_code(path)


if len(sys.argv) < 2 or len(sys.argv) > 4:
    print ("Usage: " + sys.argv[0] + " <dirname/filename> [.ext -- only file with this extension will be modified; 'all' for all files] ['only-show' -- only show filepaths that will be copyrighted]")
    sys.exit(1)

start_dir = sys.argv[1]

if not os.path.exists(start_dir):
    print("Error! " + start_dir + " does not exist!")
    sys.exit(1)   

arg = None
if len(sys.argv) >= 3:
    arg = sys.argv[2]
    if arg == "all":
        arg = None

if len(sys.argv) >= 4:
    only_show = True

if os.path.isfile(start_dir):
    visit(arg, '', [start_dir])
else:
    os.path.walk(start_dir, visit, arg)
