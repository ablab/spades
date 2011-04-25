#include "analyses.h"

void Analyses::init() {
    FILE *fp;
    char ch;
    std::string key;
    if  ((fp = fopen (m_filename,"r")) != NULL) {
        while  ((ch = getc (fp)) != EOF) {
            if (ch == '\n') {
                continue;
            }
            key += ch;
            if (key.length() == m_mer) {
                if (m_data.count(key) > 0) {
                    m_data[key] += 1;
                } else {
                    m_data[key] = 1;
                }
                key = key.substr(1);
            }
        }
        fclose (fp);
        createDatFile();
        paint();
    }
    else {
        std::cerr << "Can not open file " << m_filename;
    }
}

void Analyses::initFastaTq() {
    FILE *fp;
    bool start = false, end = false;
    char ch;
    std::string key;
    if  ((fp = fopen (m_filename,"r")) != NULL) {

        while  ((ch = getc (fp)) != EOF) {
            if (ch == '\n') {
                continue;
            }
            key += ch;
            if (key.length() == m_mer) {
                if (m_data.count(key) > 0) {
                    m_data[key] += 1;
                } else {
                    m_data[key] = 1;
                }
                key = key.substr(1);
            }
            std::cout << ch << std::endl;
        }
        fclose (fp);
        //createDatFile();
        //paint();
    }
    else {
        std::cerr << "Can not open file " << m_filename;
    }
}

void Analyses::createDatFile() {
    std::ofstream out("daily.dat");
    if(!out) {
        std::cout << "Cannot open file.\n";
        return;
    }
    for (std::map<std::string,int>::iterator it=m_data.begin() ; it != m_data.end(); ++it) {
        out << (*it).first << " " << (*it).second << "\n";
    }
    return;
}
