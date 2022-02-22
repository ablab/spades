//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace fs {
std::vector<std::filesystem::path> glob(const std::string &pattern);
}  // namespace fs
