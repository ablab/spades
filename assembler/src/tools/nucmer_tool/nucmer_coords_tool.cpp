//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * nucmer_coords_tool.cpp
 *
 *  Created on: 22.08.2011
 *
 */

#include "nucmer_coords_tool.hpp"
#include <stdio.h>
#include <cstdlib>
#include <iostream>




std::vector<NucmerAllign> ReadCoordsFile(std::string file_name){
    FILE *f = fopen(file_name.c_str(),"r");
    std::vector<NucmerAllign> res;
    if (f == NULL) {
        std::cerr<<"File not found: "<<file_name;
        return res;
    }
    char tmp_tag[1000];
    fgets(tmp_tag, 1000, f);
    while (tmp_tag[0]!='='||tmp_tag[1]!='='  ){
        fgets(tmp_tag, 1000, f);
    }
    while(true){
        NucmerAllign readedNucmerAllign;
        int count = fscanf(f, "%i %i | %i %i | %i %i | %lf |", &readedNucmerAllign.inReferenceStart, &readedNucmerAllign.inReferenceEnd
                                                 , &readedNucmerAllign.inContigStart, &readedNucmerAllign.inContigEnd
                                                 , &readedNucmerAllign.inReferenceLength, &readedNucmerAllign.inContigLength
                                                 , &readedNucmerAllign.ReferenceQuality);
        if (count != 7) break;
        fgets(tmp_tag, 1000, f);
        readedNucmerAllign.contigName = tmp_tag;
        res.push_back(readedNucmerAllign);
    }
    return res;
}
