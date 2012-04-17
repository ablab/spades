#!/usr/bin/env python

import os
import shutil
import sys

script_comment = [
    '############################################################################',
    '# Copyright (c) 2011-2012 Saint-Petersburg Academic University',
    '# All Rights Reserved',
    '# See file LICENSE for details.',
    '############################################################################', 
    '']

code_comment = [
    '//***************************************************************************',
    '//* Copyright (c) 2011-2012 Saint-Petersburg Academic University',
    '//* All Rights Reserved',
    '//* See file LICENSE for details.',
    '//****************************************************************************',
    '']

def insert_in_script(filename):
    lines = open(filename).readlines()    

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
    lines = open(filename).readlines()    

    modified = open(filename, 'w')
        
    for com_line in code_comment:
        modified.write(com_line + '\n')
    
    for line in lines:
        modified.write(line)

    modified.close()

def visit(arg, dirname, names):
    for name in names:
        ext = os.path.splitext(name)[1]
        if arg and ext != arg:
            continue
        if (ext in ['.py', '.sh']) or name.lower().startswith('cmake'):
            insert_in_script(os.path.join(dirname, name))
        elif ext in ['.hpp', '.cpp']:
            insert_in_code(os.path.join(dirname, name))

if len(sys.argv) < 2:
    print ("Usage: " + sys.argv[0] + " <src folder> [.ext -- only file with this extension will be modified]")
    sys.exit(1)

start_dir = sys.argv[1]
if not os.path.isdir(start_dir):
    print("Error! " + start_dir + " is not a directory!")
    sys.exit(1)   

arg = None
if len(sys.argv) == 3:
    arg = sys.argv[2]
os.path.walk(start_dir, visit, arg)
