#pragma once

#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <iostream>

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

void MakeDirPath(const string& path) {
  size_t slash_pos = 0;
  while ((slash_pos = path.find_first_of('/', slash_pos + 1)) != string::npos) {
    make_dir(path.substr(0, slash_pos));
  }
}

std::string GenMD5FromFiles(const std::vector<std::string> &paths) {
  std::vector<std::string> paths_s = paths;
  std::sort(paths_s.begin(), paths_s.end());

  std::string accum_string = "";
  for (auto it = paths_s.begin(); it != paths_s.end(); ++it) {
    accum_string += *it;
    accum_string += " ";
  }

  FILE *md5_output = popen(("head -n 1000 " + accum_string + " | md5sum").c_str(), "r");
  VERIFY(md5_output != NULL);

  char buf[20];
  fscanf(md5_output, "%s", buf);
  pclose(md5_output);

  return std::string(buf);
}

bool NeedToUseLongSeq(unsigned k) {
  return k > 101;
}

}

}
