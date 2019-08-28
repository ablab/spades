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
    IdType &operator[](size_t id) {
        return id_map_[id];
    }

    const IdType &operator[](size_t id) const {
        auto i = id_map_.find(id);
        VERIFY_MSG(i != id_map_.end(), "No mapped id for " << id);
        return i->second;
    }

    bool count(size_t id) const {
        return id_map_.count(id);
    }

    size_t size() const {
        return id_map_.size();
    }

private:
    std::unordered_map <size_t, IdType> id_map_;
};

//template<typename T>
//using EdgeMapper = IdMapper<typename T::EdgeId>;

}
