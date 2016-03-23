//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <iostream>

#include "compare_standard.hpp"
#include "pipeline/graphio.hpp"

namespace cap {

namespace utils {

class TmpFolderFixture {
  std::string tmp_folder_;

 public:
    TmpFolderFixture(std::string tmp_folder) {
    tmp_folder_ = tmp_folder;
    INFO("Creating " << tmp_folder_ << ": " << make_dir(tmp_folder_));
  }

  void Stub() {
    if (1 + 1 == 1) {
      INFO("LOL WAT");
    }
  }

  ~TmpFolderFixture() {
    // Deleting temporary dir iff all indices are destructed
    int ret_code = rmdir(tmp_folder_.c_str());
    INFO("Removing temporary : " << (ret_code ? "failure, some indices were \
          not deleted" : "success"));
  }
};

int add_time(double &time, int multiplier = 1, int ret = 0) {
  timeval curtime;
  gettimeofday(&curtime, 0);
  time += multiplier * (double(curtime.tv_sec) + double(curtime.tv_usec) * 1e-6);
  return ret;
}

//todo remove
void MakeDirPath(const std::string& path) {
  if (path.size() == 0) {
    TRACE("Somewhat delirium: trying to create directory ``");
    return;
  }

  size_t slash_pos = 0;
  while ((slash_pos = path.find_first_of('/', slash_pos + 1)) != std::string::npos) {
    make_dir(path.substr(0, slash_pos));
  }
  if (path[path.size() - 1] != '/') {
    make_dir(path);
  }
}

//todo remove
bool DirExist(std::string path) {
  struct stat st;
  return (stat(path.c_str(), &st) == 0) && (S_ISDIR(st.st_mode));
}

ContigStreams OpenStreams(const vector<string>& filenames) {
  ContigStreams streams;
  for (auto it = filenames.begin(); it != filenames.end(); ++it) {
    DEBUG("Opening stream from " << *it);
    streams.push_back(make_shared<io::FileReadStream>(*it));
  }
  return streams;
}

std::string GetMD5CommandString() {
  static std::string answer = "";

  if (answer != "") {
    return answer;
  }

  FILE *output;
  char buf[40];
  output = popen("echo a | md5sum 2> /dev/null", "r");
  if (1 == fscanf(output, "%s", buf) &&
      strcmp(buf, "60b725f10c9c85c70d97880dfe8191b3") == 0) {
    return answer = "md5sum ";
  }
  pclose(output);

  output = popen("echo a | md5 2> /dev/null", "r");
  if (1 == fscanf(output, "%s", buf) &&
      strcmp(buf, "60b725f10c9c85c70d97880dfe8191b3") == 0) {
    return answer = "md5 ";
  }
  pclose(output);

  return answer = "head -c 20 ";

}

std::string GenMD5FromFiles(const std::vector<std::string> &paths,
                            const std::string &salt = "") {
  VERIFY(!paths.empty());
  std::vector<std::string> paths_s = paths;
  //std::sort(paths_s.begin(), paths_s.end());

  std::string accum_string = "";
  for (auto it = paths_s.begin(); it != paths_s.end(); ++it) {
    accum_string += *it;
    accum_string += " ";
  }

  std::cerr << "Using " << GetMD5CommandString() << std::endl;

  FILE *md5_output = popen(("(head -n 1000 " + accum_string + "&& echo " + salt
                                + ") | " + GetMD5CommandString()).c_str(), "r");
  VERIFY(md5_output != NULL);

  char buf[40];
  VERIFY(1 == fscanf(md5_output, "%s", buf));
  pclose(md5_output);

  return std::string(buf);
}

bool NeedToUseLongSeq(size_t k) {
  return k > 99;
}

}

}
