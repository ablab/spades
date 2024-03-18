//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
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
