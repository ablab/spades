#ifndef ANALYSES_H
#define ANALYSES_H

#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <fstream> 
#include <map>

class Analyses {
private:
    std::map <std::string, int> m_data;
    int m_mer;
    char *m_filename;

    void createDatFile();

    void paint() {
        system("gnuplot ./file.gnu");
    }
    void init();
public:
    Analyses(char *argv, char *mer) : m_filename(argv), m_mer(atoi(mer)) {
        init();
    }
};

#endif
