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

namespace binary {

class BinSaveFile {
public:
    BinSaveFile(const std::string &filename)
            : str_(filename, std::ios::binary) {
    }

    template<typename T>
    BinSaveFile &operator<<(const T &value) {
        io::binary::BinWrite(str_, value);
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

class BinLoadFile {
public:
    BinLoadFile(const std::string &filename)
            : str_(filename, std::ios::binary) {
    }

    template<typename T>
    BinLoadFile &operator>>(T &value) {
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

    std::istream &stream() {
        return str_;
    }

private:
    std::ifstream str_;
};

/**
 * @brief  An interface that can consistently save and load some component T.
 * @tparam Env (optional) parameters to be used in loading (typically an EdgeId mapper).
 */
template<typename T, typename... Env>
class IOBase {
    virtual void Save(const std::string &basename, const T &value) = 0;
    virtual bool Load(const std::string &basename, T &value, const Env &... env) = 0;
};

/**
  * @brief  IOTraits<T>::type is the concrete class is needed for (de)serialization of some type T.
  */
template<typename T>
struct IOTraits;

/**
 * @brief  An abstract saver/loader which uses a single file for its component.
 */
template<typename T, typename... Env>
class IOSingle : public IOBase<T, Env...> {
public:


    void Save(const std::string &basename, const T &value) override {
        std::string filename = basename + this->ext_;
        BinSaveFile file(filename);
        DEBUG("Saving " << this->name_ << " into " << filename);
        VERIFY(file);
        this->SaveImpl(file, value);
    }

    bool Load(const std::string &basename, T &value, const Env &... env) override {
        std::string filename = basename + this->ext_;
        BinLoadFile file(filename);
        DEBUG("Loading " << this->name_ << " from " << filename);
        if (!file)
            return false;
        this->LoadImpl(file, value, env...);
        return true;
    }

protected:
    const char *name_, *ext_;

    IOSingle(const char *name, const char *ext)
            : name_(name), ext_(ext) {
    }

    virtual void SaveImpl(BinSaveFile &file, const T &value) = 0;
    virtual void LoadImpl(BinLoadFile &file, T &value, const Env &... env) = 0;

private:
    DECL_LOGGER("BinaryIO");
};

/**
 * @brief  A default implementor of single-file binary saving/loading.
 */
template<typename T, typename... Env>
class IOSingleDefault : public IOSingle<T, Env...> {
public:
    IOSingleDefault(const char *name, const char *ext)
            : IOSingle<T, Env...>(name, ext) {
    }

protected:
    void SaveImpl(BinSaveFile &file, const T &value) override {
        file << value;
    }

    void LoadImpl(BinLoadFile &file, T &value, const Env &... env) override {
        value.BinRead(file.stream(), env...);
    }
};

/**
 * @brief  An saver/loader over a vector-like collection (one file with numerical suffix per element).
 */
template<typename T, typename Env>
class IOCollection : public IOBase<T, Env> {
public:
    typedef std::unique_ptr<IOSingle<typename T::value_type, Env>> IOPtr;

    IOCollection(IOPtr io)
        : io_(std::move(io)) {
    }

    void Save(const std::string &basename, const T &value) override {
        for (size_t i = 0; i < value.size(); ++i) {
            io_->Save(basename + "_" + std::to_string(i), value[i]);
        }
    }

    bool Load(const std::string &basename, T &value, const Env &env) override {
        bool res = true;
        for (size_t i = 0; i < value.size(); ++i) {
            res &= io_->Load(basename + "_" + std::to_string(i), value[i], env);
        }
        return res;
    }

private:
    IOPtr io_;
};

} // namespace binary

} //namespace io
