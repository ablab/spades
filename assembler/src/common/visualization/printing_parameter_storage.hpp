#pragma once

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/components/graph_component.hpp"
#include <boost/optional.hpp>
#include <map>

using namespace omnigraph;

namespace visualization {

namespace printing_parameter_storage {

template<typename ElementId, typename Value>
class ParameterStorage {
public:
    virtual Value GetValue(ElementId element) const = 0;

    virtual ~ParameterStorage() {
    }
};

template<typename ElementId, typename Value>
class MapParameterStorage : public virtual ParameterStorage<ElementId, Value> {
public:
    typedef std::map<ElementId, Value> Storage;
private:
    template<class It>
    static std::map<ElementId, std::string> ConstructMap(It begin, It end, const std::string &color) {
        std::map<ElementId, std::string> result;
        for (auto it = begin; it != end; ++it) {
            result.insert({*it, color});
        }
        return result;
    }

protected:
    std::map<ElementId, Value> storage_;
private:
    boost::optional<Value> default_value_;
public:
    MapParameterStorage(const std::string &default_value):
            default_value_(default_value) {
    }

    MapParameterStorage(const Storage &storage, Value default_value = Value()):
            storage_(storage), default_value_(std::move(default_value)) {
    }

    template<class It>
    MapParameterStorage(It begin, It end, const Value &value, const std::string &default_value):
            storage_(ConstructMap(begin, end, value)), default_value_(default_value) {
    }


    Value GetValue(ElementId element) const {
        auto it = storage_.find(element);
        if (it == storage_.end()) {
            VERIFY(default_value_);
            return default_value_.get();
        }
        return it->second;
    }
};

template<typename ElementId, typename Value>
class DecoratorParameterStorage : public virtual ParameterStorage<ElementId, Value> {
private:
    ParameterStorage<ElementId, Value> inner_storage_;
public:
    DecoratorParameterStorage(ParameterStorage<ElementId, Value> inner_storage) : inner_storage_(
            inner_storage) {
    }

    Value GetInnerValue(ElementId element) {
        return inner_storage_.GetValue(element);
    }
};

}
}