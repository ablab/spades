//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <exception>
#include <iostream>
int pathracer_main(int argc, char* argv[]);

int main(int argc, char* argv[]) {
    try {
        return pathracer_main(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown object caught" << std::endl;
        return 1;
    }
}
