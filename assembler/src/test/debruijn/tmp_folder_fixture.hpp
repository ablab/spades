//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/filesystem/path_helper.hpp"
#include <string>

class TmpFolderFixture {
    std::string tmp_folder_;

public:
    TmpFolderFixture(std::string tmp_folder = "tmp")
            : tmp_folder_(tmp_folder) {
        fs::make_dirs(tmp_folder_);
    }

    ~TmpFolderFixture() {
        fs::remove_dir(tmp_folder_);
    }

    const std::string &tmp_folder() const {
        return tmp_folder_;
    }    
};
