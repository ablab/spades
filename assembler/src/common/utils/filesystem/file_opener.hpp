//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once 

#include <fstream>
#include <string>

/// @returns an opened file stream
/// @throw std::ios_base::failure if the file does not exists
/// @note Be careful with std::ios_base::failbit and reading a file until the eof
template <class FileName>
std::ifstream open_file(FileName && fileName,
                       std::ios_base::openmode mode = std::ios_base::in,
                       std::ios_base::iostate exceptionBits = std::ios_base::failbit | std::ios_base::badbit)
{
    std::ifstream file(std::forward<FileName>(fileName), mode);
    if (!file.is_open())
        throw std::ios_base::failure("Cannot open file '" + std::string(fileName) + '\'');
    file.exceptions(exceptionBits);
    return file;
}