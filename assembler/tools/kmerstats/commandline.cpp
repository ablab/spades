#include "commandline.h"

#include <cstdlib>
#include <getopt.h>
#include <iostream>

#include "analyses.h"

CommandLine::CommandLine(char **argv, int argc) {
    int c;

    while (1) {
        int option_index = 0;
        static struct option long_options[] = {
            {"file", 1, 0, 0},
            {"kmer", 1, 0, 0},
            {"help", 0, 0, 0},
            {"merge", 0, 0, 0},
            {0, 0, 0, 0}
        };

        c = getopt_long (argc, argv, "f:k:h",
                         long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
        case 0:
            if (long_options[option_index].name == "file") {
                m_line = optarg;
            } else if (long_options[option_index].name == "kmer") {
                m_kmer = optarg;
            }
            break;

        case 'f':
            m_line = optarg;
            break;

        case 'k':
            m_kmer = optarg;
            break;

        case 'h':
        default:
            std::cout << "Help: ./kmerstat -f <filename> -k <kmer>" << std::endl;
        }
    }

    Analyses a(m_line, m_kmer);
}

