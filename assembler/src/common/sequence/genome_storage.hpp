//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include "sequence.hpp"
#include "nucl.hpp"

struct Chromosome{
    std::string name;
    std::string sequence;
    Chromosome(string chr_name, string seq): name(chr_name), sequence(seq){}
};

class GenomeStorage {
//all chromosomes glued
    std::string s_;
    std::vector<Chromosome> full_genome_;

    std::string ACGTise(const std::string &s) const {
        std::stringstream ss;
        for (size_t i = 0; i < s.length(); i++){
            if (!is_nucl(s[i])) {
                char upper_cased = char(std::toupper(s[i]));
                if (is_nucl(upper_cased))
                    ss << upper_cased;
            } else
                ss << s[i];
        }
        return ss.str();
    }
public:
    GenomeStorage() {
    }

    GenomeStorage(const std::string &s): s_(s), full_genome_() {
        full_genome_.push_back(Chromosome("genome", ACGTise(s_)));
    }

    GenomeStorage(const vector<Chromosome> &chromosomes): full_genome_(chromosomes) {
        std::stringstream ss;
        for (const auto &s: chromosomes) {
            ss << ACGTise(s.sequence);
//do we need a separator between?
        }
        s_ = ss.str();
    }

    GenomeStorage(const vector<string> &chromosomes): full_genome_() {
        std::stringstream ss;
        int count = 0;
        for (const auto &s: chromosomes) {
            count ++;
            std::string fxd = ACGTise(s);
            full_genome_.push_back(Chromosome("chr" + std::to_string(count), fxd));
            ss << fxd;
//do we need a separator between?
        }
        s_ = ss.str();
    }


    //TODO exterminate this where possible
    Sequence GetSequence() const {
        stringstream ss;
        size_t l = 0, r = 0;
        for(size_t i = 0; i < s_.size(); i++) {
            if (!is_nucl(s_[i]) ) {
                if (r > l) {
                    ss << s_.substr(l, r - l);
                }
                r = i + 1;
                l = i + 1;
            } else {
                r++;
            }
        }
        if (r > l) {
            ss << s_.substr(l, r - l);
        }
        return Sequence(ss.str());
    }

    std::vector<Chromosome> GetChromosomes() const{
        return full_genome_;
    }

    void SetSequence(const Sequence &s) {

        s_ = s.str();
    }

    std::string str() const {
        return s_;
    }

    size_t size() const {
        return s_.size();
    }
};

