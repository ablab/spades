//
// Created by lab42 on 8/19/15.
//

#ifndef GENOME_STORAGE_HPP_
#define GENOME_STORAGE_HPP_

#include <string>
#include "sequence/sequence.hpp"
namespace debruijn_graph {
    class GenomeStorage {
    private:
        std::string s_;
    public:
        GenomeStorage():s_(""){
        }

        GenomeStorage(const std::string &s): s_(s){
        }

        Sequence GetSequence() const;
        size_t size() const;
    };
}
#endif //PROJECT_GENOME_STORAGE_HPP
