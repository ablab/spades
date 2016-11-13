//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * qual.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef QUAL_HPP_
#define QUAL_HPP_

#include <string>
//todo really strange class
class Quality {
public:

    Quality(const std::string &s) : qual_(s) {
    }

    int operator[](size_t i) const {
        return qual_[i];
    }

    std::string str() const { // copying (defensive)!
        return qual_;
    }

private:
    std::string qual_;
  //friend class ireadstream;
};

#endif /* QUAL_HPP_ */
