//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/verify.hpp"

#include <unordered_map>

namespace io {

template<typename IdType>
class IdMapper {
public:
    void map(const IdType &from, size_t to) {
        id_map_.emplace(std::make_pair(to, from));
        back_map_.emplace(std::make_pair(from, to));
    }

    const IdType &operator[](size_t id) const {
        auto i = id_map_.find(id);
        VERIFY_MSG(i != id_map_.end(), "No mapped id for " << id);
        return i->second;
    }

    size_t operator[](const IdType &id) const {
        auto i = back_map_.find(id);
        VERIFY_MSG(i != back_map_.end(), "No mapped id for " << id);
        return i->second;
    }

    bool count(size_t id) const {
        return id_map_.count(id);
    }

    size_t size() const {
        return id_map_.size();
    }

private:
    std::unordered_map<size_t, IdType> id_map_;
    std::unordered_map<IdType, size_t> back_map_;
};

//template<typename T>
//using EdgeMapper = IdMapper<typename T::EdgeId>;

}
