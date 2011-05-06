#ifndef ANALYSES_H
#define ANALYSES_H

#include <iostream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <fstream> 
#include <map>
#include "algoritms/SpaceSaving.h"
#include "algoritms/CountMinSketch.h"
#include "algoritms/Frequent.h"
#include "algoritms/LossyCounting.h"
#include "algoritms/DistinctElements.h"

class Analyses {
private:
    std::map <std::string, int> m_data;
    int m_mer;
    char *m_filename;
    char *m_algname;

    std::map <std::string,  void (*)(char *, off_t, int, double, double) > m_algorithm;
    void createDatFile();

    void paint() {
        system("gnuplot ./file.gnu");
    }
    void init();
    void initFastTq();
public:
    Analyses(char *argv, char *mer, char *alg = "") : m_filename(argv), m_mer(atoi(mer)), m_algname(alg) {
        m_algorithm["spacesave"] = WrapperSpaceSaving;
        m_algorithm["cms"] = WrapperCountMinSketch;
        m_algorithm["frequent"] = WrapperFrequent;
        m_algorithm["lossy"] = WrapperLossyCounting;
        m_algorithm["distinct"] = WrapperDistinctElements;

        if ((strstr(argv, ".fastatq")) != NULL) {
            initFastTq();
        } else {
            init();
        }
    }
};

#endif
