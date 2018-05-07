//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace io {

class SaveFile {
public:
    SaveFile(const std::string &filename)
            : str_(filename, std::ios::binary) {
    }

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
    LoadFile(const std::string &filename)
            : str_(filename, std::ios::binary) {
    }

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

template<typename T>
class IOBase {
public:
    virtual void Save(const std::string &basename, const T &value) = 0;
    virtual void Load(const std::string &basename, T &value) = 0;
};

template <typename T>
class IOSingle : public IOBase<T> {
public:
    void Save(const std::string &basename, const T &value) override {
        std::string filename = basename + ext_;
        SaveFile file(filename);
        DEBUG("Saving " << name_ << " into " << filename);
        VERIFY(file);
        this->SaveImpl(file, value);
    }

    void Load(const std::string &basename, T &value) override {
        std::string filename = basename + ext_;
        LoadFile file(filename);
        DEBUG("Loading " << name_ << " from " << filename);
        VERIFY(file);
        this->LoadImpl(file, value);
    }

protected:
    IOSingle(const char *name, const char *ext)
            : name_(name), ext_(ext), version_(0) {
    }
    const char *name_, *ext_;
    unsigned version_;

private:
    virtual void SaveImpl(SaveFile &file, const T &value) = 0;
    virtual void LoadImpl(LoadFile &file, T &value) = 0;
};

template<typename T>
class IOCollection : public IOBase<T> {
public:
    typedef std::unique_ptr<IOBase<typename T::value_type>> IOPtr;

    IOCollection(IOPtr io)
        : io_(std::move(io)) {
    }

    void Save(const std::string &basename, const T &value) override {
        for (size_t i = 0; i < value.size(); ++i) {
            io_->Save(basename + "_" + std::to_string(i), value[i]);
        }
    }

    void Load(const std::string &basename, T &value) override {
        for (size_t i = 0; i < value.size(); ++i) {
            io_->Load(basename + "_" + std::to_string(i), value[i]);
        }
    }

private:
    IOPtr io_;
};

}
