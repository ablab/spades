//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <iostream>
#include <unordered_map>
#include <cstring>
#include <string>
#include <vector>

#include "sequence/sequence.hpp"

namespace cap {

class Serializer {
  public:
    Serializer(std::ostream &os) : os_(os) {
    }

    template <class T>
    void WriteLine(const std::string &key, const T &value) {
        os_ << key << kDelimiter;
        SerializeToStream(value, os_);
        os_ << std::endl;
    }

    template <class T>
    std::string Serialize(const T &value) const {
        std::stringstream ss;
        SerializeToStream(value, ss);
        return ss.str();
    }

  private:
    static const char kDelimiter = char(1);

    template <class T>
    void SerializeToStream(const T &value, std::ostream &out) const {
        out << value;
        out << kDelimiter;
    }
    template <class T>
    void SerializeToStream(const std::vector<T> &value, std::ostream &out) const;
    template<class T1, class T2>
    void SerializeToStream (const std::pair<T1, T2> &value, std::ostream &out) const;

    std::ostream &os_;
};

template<class T>
void Serializer::SerializeToStream(const std::vector<T> &value, std::ostream &out) const {
    //out << value.size() << kDelimiter;
    SerializeToStream(value.size(), out);
    for (const auto &x : value) {
        SerializeToStream(x, out);
    }
}

template<class T1, class T2>
void Serializer::SerializeToStream (const std::pair<T1, T2> &value, std::ostream &out) const {
    SerializeToStream(value.first, out);
    SerializeToStream(value.second, out);
}

class Deserializer {
  public:
    Deserializer(std::istream &is) : is_(is) {
    }

    void ReadStream() {
        while (!is_.eof()) {
            std::string key;
            std::string value;

            std::getline(is_, key, kDelimiter);
            std::getline(is_, value);
            read_map_[key] = value;

        }
    }

    template <class T>
    void ReadValue(const std::string &key, T &value) const {
        if (read_map_.count(key) == 0)
            return;
        Deserialize(read_map_.at(key), value);
    }

    template <class T>
    void Deserialize(const std::string &s, T &value) const {
        std::stringstream ss(s);
        DeserializeFromStream(ss, value);
    }

  private:
    static const char kDelimiter = char(1);

    template <class T>
    void DeserializeFromStream(std::istream &is, T &value) const {
        is >> value;
        is.ignore(1, kDelimiter);
    }
    template <class T>
    void DeserializeFromStream(std::istream &is, std::vector<T> &value) const;
    template<class T1, class T2>
    void DeserializeFromStream (std::istream &is,
            std::pair<T1, T2> &value) const;
    void DeserializeFromStream(std::istream &is, Sequence &s) const;

    std::istream &is_;
    std::unordered_map<std::string, std::string> read_map_;
};

template <>
void Deserializer::DeserializeFromStream<std::string>(std::istream &is,
        std::string &value) const {
    std::getline(is, value, kDelimiter);
}

template <class T>
void Deserializer::DeserializeFromStream(std::istream &is,
        std::vector<T> &value) const {
    size_t size = 0;
    DeserializeFromStream(is, size);
    value.resize(size);

    for (size_t i = 0; i < size; ++i) {
        DeserializeFromStream(is, value[i]);
    }
}

template<class T1, class T2>
void Deserializer::DeserializeFromStream (std::istream &is,
        std::pair<T1, T2> &value) const {
    DeserializeFromStream(is, value.first);
    DeserializeFromStream(is, value.second);
}

void Deserializer::DeserializeFromStream(std::istream &is,
        Sequence &s) const {
    std::string str;
    DeserializeFromStream(is, str);
    s = Sequence(str.c_str());
}

}
