#include "analyses.h"
#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

void Analyses::init() {
    std::string key;

    struct stat sb;
    off_t len;
    char *p;
    int fd;

    fd = open (m_filename, O_RDONLY);
    if (fd == -1) {
        perror ("Error: could not open");
        return;
    }

    if (fstat (fd, &sb) == -1) {
        perror ("Error: have a fstat problem");
        return;
    }

    if (!S_ISREG (sb.st_mode)) {
        fprintf (stderr, "%s is not a file\n", m_filename);
        return;
    }

    p = (char *)mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (p == MAP_FAILED) {
        perror ("mmap");
        return;
    }

    if (close (fd) == -1) {
        perror ("close");
        return;
    }

    for (len = 0; len < sb.st_size; len++) {
        if (p[len] == '\n') {
            continue;
        }
        key += p[len];
        if (key.length() == m_mer) {
            if (m_data.count(key) > 0) {
                m_data[key] += 1;
            } else {
                m_data[key] = 1;
            }
            key = key.substr(1);
        }
    }

    if (munmap (p, sb.st_size) == -1) {
        perror ("munmap");
        return;
    }
    createDatFile();
    paint();
}

void Analyses::initFastTq() {
    FILE *fp;
    char ch;
    char buffer[100];	
    std::string key;

    if  ((fp = fopen (m_filename,"r")) != NULL) {
	int count = 1;	
	while  ( ! feof (fp)) {
	    fgets (buffer, 100, fp);
	    if (count == 2) {
    		int i = 0;
		while (buffer[i] != '\n') {
                    key += buffer[i];
                    if (key.length() == m_mer) {
                        if (m_data.count(key) > 0) {
                            m_data[key] += 1;
                        } else {
                            m_data[key] = 1;
                        }
                        key = key.substr(1);
                    }
                    ++i;
		}
            } 
	    if (count <= 3) { 	
	    	++count; 	
	    } else {
		count = 1;
	    }		
        }
        fclose (fp);

	for (std::map<std::string,int>::iterator it=m_data.begin() ; it != m_data.end(); it++ )
	    std::cout << (*it).first << " => " << (*it).second << std::endl;
        createDatFile();
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
