#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#include <string>
#include <iostream>

class CommandLine {
    char *m_line;
    char *m_kmer;
public:
    CommandLine(char **line, int count);
};


#endif // COMMANDLINE_H
