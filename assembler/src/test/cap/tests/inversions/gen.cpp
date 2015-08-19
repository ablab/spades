//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <cstdio>
#include <cstring>
#include <string>
#include <algorithm>
#include <cstdlib>

const std::string atgc = "ATGC";
const double MAX_REVERSE_LENGTH = 4./5;
const size_t MIN_REVERSE_LENGTH = 30;

int int_rand() {
  return abs(rand() * RAND_MAX + rand());
}

char complementary(char c) {
  return atgc[atgc.find(c) ^ 1];
}

void introduce_inversion(std::string &s, int &pos, int &len) {
  size_t L = s.size();

  len = MIN_REVERSE_LENGTH +
      (int_rand() % (int)(MAX_REVERSE_LENGTH * L - MIN_REVERSE_LENGTH));
  pos = MIN_REVERSE_LENGTH + (int_rand() % (L - len - 2 * MIN_REVERSE_LENGTH));

  std::reverse(s.begin() + pos, s.begin() + pos + len);
  for (int i = pos; i < pos + len; ++i) {
    s[i] = complementary(s[i]);
  }
}

void write_fasta(const char *filename, const char *desc, const char *seq) {
  FILE *fd = fopen(filename, "w");
  fprintf(fd, ">%s\n", desc);
  for (size_t i = 0;; ++i) {
    if (seq[i] == 0) break;
    fputc(seq[i], fd);
    if (i % 80 == 79) fputc('\n', fd);
  }
  fputs("\n", fd);
  fclose(fd);
}

// Takes L and N --- length of the sequence and number of inversions
int main() {
  srand(2391);

  unsigned int L, N;
  char file1[30],
       file2[30];
  scanf("%d %d %s %s", &L, &N, file1, file2);

  std::pair<int, int> stats[N];
  std::string s;
  for (size_t i = 0; i < L; ++i) {
    s += atgc[rand() % 4];
  }
  std::string t = s;

  for (size_t i = 0; i < N; ++i) {
    introduce_inversion(t, stats[i].second, stats[i].first);
  }
  std::sort(stats, stats + N);
  std::reverse(stats, stats + N);

  printf("%d %d\n", L, N);
  for (size_t i = 0; i < N; ++i) {
    printf("%d %d\n", stats[i].first, stats[i].second);
  }

  write_fasta(file1, "forward", s.c_str());
  write_fasta(file2, "reversed", t.c_str());
  
  return 0;
}
