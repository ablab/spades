//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/binary.hpp"
#include "utils/logger/logger.hpp"
#include "utils/filesystem/path_helper.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace io {

namespace binary {

/**
 * @brief  An interface that can consistently save and load some component T.
 */
template<typename T>
struct IOBase {
    virtual void Save(const std::string &basename, const T &value) = 0;
    virtual bool Load(const std::string &basename, T &value) = 0;
    virtual ~IOBase() {}
};

/**
 * @brief  IOTraits<T>::type is the concrete class is needed for (de)serialization of some type T.
 */
template<typename T>
struct IOTraits;

/**
 * @brief  A convenient template function that saves a component into a file
 *         calling an appropriate ComponentIO.
 */
template<typename T>
void Save(const std::string &basename, const T &value) {
    typename IOTraits<T>::Type io;
    io.Save(basename, value);
}

/**
 * @brief  A convenient template function that reads a component from a file
 *         calling an appropriate ComponentIO.
 */
template<typename T>
bool Load(const std::string &basename, T &value) {
    typename IOTraits<T>::Type io;
    return io.Load(basename, value);
}

/**
 * @brief  An interface that can consistently (de)serialize some component T using a binary stream wrapper.
 */
template<typename T>
struct IOStream : public IOBase<T> {
    virtual void Write(BinOStream &stream, const T &value) = 0;
    virtual void Read(BinIStream &stream, T &value) = 0;
};

/**
 * @brief  A convenient template function that serializes a component into an STL ostream
 *         calling an appropriate ComponentIO.
 */
template<typename T>
void Write(std::ostream &str, const T &value) {
    BinOStream str_(str);
    typename IOTraits<T>::Type io;
    io.Write(str_, value);
}

/**
 * @brief  A convenient template function that deserializes a component from an STL stream
 *         calling an appropriate ComponentIO.
 */
template<typename T, typename... Env>
bool Load(std::istream &str, T &value, Env... env) {
    BinIStream str_(str);
    typename IOTraits<T>::Type io;
    return io.Read(str_, value, env...);
}

/**
 * @brief  An abstract saver/loader which uses a single file for its component.
 */
template<typename T>
class IOSingle : public IOStream<T> {
public:
    IOSingle(const char *name, const char *ext)
            : name_(name), ext_(ext) {
    }

    void Save(const std::string &basename, const T &value) override {
        std::string filename = basename + this->ext_;
        std::ofstream file(filename, std::ios::binary);
        DEBUG("Saving " << this->name_ << " into " << filename);
        VERIFY(file);
        BinOStream writer(file);
        this->Write(writer, value);
    }

    /**
     * @return false if the file is missing. true if the component was successfully loaded.
     *         Fails if the file is present but cannot be read.
     */
    bool Load(const std::string &basename, T &value) override {
        std::string filename = basename + this->ext_;
        if (!fs::check_existence(filename))
            return false;
        std::ifstream file(filename, std::ios::binary);
        VERIFY_MSG(file, "Failed to read " << filename);
        DEBUG("Loading " << this->name_ << " from " << filename);
        BinIStream reader(file);
        this->Read(reader, value);
        return true;
    }

private:
    const char *name_, *ext_;

    DECL_LOGGER("BinaryIO");
};

/**
 * @brief  A default implementor of single-file binary saving/loading.
 */
template<typename T>
class IOSingleDefault : public IOSingle<T> {
public:
    IOSingleDefault(const char *name, const char *ext)
            : IOSingle<T>(name, ext) {
    }

    void Write(BinOStream &str, const T &value) override {
        str << value;
    }

    void Read(BinIStream &str, T &value) override {
        value.BinRead(str.stream());
    }
};

/**
 * @brief  An saver/loader over a vector-like collection (one file with numerical suffix per element).
 */
template<typename T>
class IOCollection : public IOBase<T> {
public:
    typedef std::unique_ptr<IOSingle<typename T::value_type>> IOPtr;

    IOCollection(IOPtr io)
        : io_(std::move(io)) {
    }

    void Save(const std::string &basename, const T &value) override {
        for (size_t i = 0; i < value.size(); ++i) {
            io_->Save(basename + "_" + std::to_string(i), value[i]);
        }
    }

    bool Load(const std::string &basename, T &value) override {
        bool res = true;
        for (size_t i = 0; i < value.size(); ++i) {
            res &= io_->Load(basename + "_" + std::to_string(i), value[i]);
        }
        return res;
    }

private:
    IOPtr io_;
};

} // namespace binary

} //namespace io
