/*
 * seq.cpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#include "seq.hpp"

MatePair::MatePair(const std::string &s1, const std::string &s2, const int id_) : id(id_), seq1(s1), seq2(s2) {
}

MatePair::MatePair(const MatePair &mp) : id(mp.id), seq1(mp.seq1), seq2(mp.seq2) {
}
