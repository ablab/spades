//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#include <string>
#include <iostream>

class CommandLine {
    char *m_line;
    char *m_kmer;
    char *m_nameAlg;
    bool isAlg;
public:
    CommandLine(char **line, int count);
};


#endif // COMMANDLINE_H
