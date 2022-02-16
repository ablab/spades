//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>

class TmpFolderFixture {
    std::filesystem::path tmp_folder_;

public:
    TmpFolderFixture(std::filesystem::path tmp_folder = "tmp")
            : tmp_folder_(tmp_folder) {
        create_directories(tmp_folder_);
    }

    ~TmpFolderFixture() {
        remove(tmp_folder_);
    }

    const std::filesystem::path &tmp_folder() const {
        return tmp_folder_;
    }    
};
