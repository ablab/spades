//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <iostream>

#include "compare_standard.hpp"

namespace cap {
typedef io::SingleRead Contig;
typedef io::IReader<Contig> ContigStream;
typedef	io::MultifileReader<io::SingleRead> CompositeContigStream;
typedef	io::RCReaderWrapper<io::SingleRead> RCWrapper;
}

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
  time += multiplier * (curtime.tv_sec + curtime.tv_usec * 1e-6);
  return ret;
}

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

bool DirExist(std::string path) {
  struct stat st;
  return (stat(path.c_str(), &st) == 0) && (S_ISDIR(st.st_mode));
}

vector<cap::ContigStream*> OpenStreams(const vector<string>& filenames) {
  vector<ContigStream*> streams;
  for (auto it = filenames.begin(); it != filenames.end(); ++it) {
    DEBUG("Opening stream from " << *it);
    streams.push_back(new io::Reader(*it));
  }
  return streams;
}

std::string GenMD5FromFiles(const std::vector<std::string> &paths, const std::string &salt = "") {
  std::vector<std::string> paths_s = paths;
  //std::sort(paths_s.begin(), paths_s.end());

  std::string accum_string = "";
  for (auto it = paths_s.begin(); it != paths_s.end(); ++it) {
    accum_string += *it;
    accum_string += " ";
  }

  FILE *md5_output = popen(("(head -n 1000 " + accum_string + "&& echo " + salt + ") | md5sum").c_str(), "r");
  VERIFY(md5_output != NULL);

  char buf[20];
  fscanf(md5_output, "%s", buf);
  pclose(md5_output);

  return std::string(buf);
}

bool NeedToUseLongSeq(unsigned k) {
  return k > 99;
}

}

}
