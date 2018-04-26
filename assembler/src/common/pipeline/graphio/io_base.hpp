//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace debruijn_graph {

namespace graphio {

class SaveFile {
public:
    SaveFile(const std::string &filename) :
            str_(filename, std::ios::binary) {}

    template<typename T>
    SaveFile &operator<<(const T &value) {
        io::binary::BinWrite(str_, value);
        //VERIFY(!str.fail());
        return *this;
    }

    operator bool() const {
        return (bool) str_;
    }

    bool operator!() const {
        return !str_;
    }

private:
    std::ofstream str_;
};

class LoadFile {
public:
    LoadFile(const std::string &filename) :
            str_(filename, std::ios::binary) {}

    template<typename T>
    LoadFile &operator>>(T &value) {
        io::binary::BinRead(str_, value);
        //VERIFY(!str.fail());
        return *this;
    }

    template<typename T>
    T Read() {
        T result;
        (*this) >> result;
        return result;
    }

    operator bool() const {
        return (bool) str_;
    }

    bool operator!() const {
        return !str_;
    }

    std::istream &str() { //TODO: get rid of this hack
        return str_;
    }

private:
    std::ifstream str_;
};

template <typename T>
class IOBase {
public:
    typedef T Type;

    void Save(const std::string &basename, const T &value) {
        std::string filename = basename + ext_;
        SaveFile file(filename);
        DEBUG("Saving " << name_ << " into " << filename);
        VERIFY(file);
        this->SaveImpl(file, value);
    }

    void Load(const std::string &basename, T &value) {
        std::string filename = basename + ext_;
        LoadFile file(filename);
        DEBUG("Loading " << name_ << " from " << filename);
        VERIFY(file);
        this->LoadImpl(file, value);
    }

protected:
    IOBase(const char *name, const char *ext):
            name_(name), ext_(ext), version_(0) {}
private:
    virtual void SaveImpl(SaveFile &file, const T &value) = 0;
    virtual void LoadImpl(LoadFile &file, T &value) = 0;

    const char *name_, *ext_;
    unsigned version_;
};

}

}
