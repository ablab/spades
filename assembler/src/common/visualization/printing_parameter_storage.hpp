#pragma once

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "utils/standard_base.hpp"
#include "assembly_graph/components/graph_component.hpp"

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
private:
private:
    template<class It>
    static map<ElementId, string> ConstructMap(It begin, It end, const string &color) {
        map<ElementId, string> result;
        for (auto it = begin; it != end; ++it) {
            result.insert(make_pair(*it, color));
        }
        return result;
    }

protected:
    map<ElementId, Value> storage_;
private:
    boost::optional<Value> default_value_;
public:
    MapParameterStorage(const string &default_value) : default_value_(default_value) {
    }

    MapParameterStorage(map<ElementId, Value> storage, Value default_value) : storage_(storage),
                                                                              default_value_(default_value) {
    }

    MapParameterStorage(map<ElementId, Value> storage) : storage_(storage) {
    }

    template<class It>
    MapParameterStorage(It begin, It end, const Value &value, const string &default_value) : storage_(
            ConstructMap(begin, end, value)), default_value_(default_value) {
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