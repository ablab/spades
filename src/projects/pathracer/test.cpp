#include <gtest/gtest.h>

#include "find_best_path.hpp"
#include "fees.hpp"

double levenshtein_string_score(const std::string &s, const std::string &query) {
  auto fees = hmm::levenshtein_fees(query);
  fees.minimal_match_length = 0;
  return -score_sequence(fees, s);
}

double levenshtein_substring_score(const std::string &s, const std::string &query) {
  auto fees = hmm::levenshtein_fees(query);
  fees.minimal_match_length = 0; // TODO Automatically reduce required mathed string length for short HMM
  return -score_subsequence(fees, s);
}

TEST(LevenshteinZero, LEVENSHTEIN_SUBSTRING) {
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("A", "A"), 0);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AAAAA", "A"), 0);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AAAAACGTAAAAAAACGT", "AAA"), 0);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AAAAACGTAAAAAAACGT", ""), 0);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("", ""), 0);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AAAAACGTAAAAAAACGT", "CGT"), 0);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AAAAACGTAAAAAAACGT", "AAAAACGTAAAAAAACGT"), 0);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("TTTTTTTTTTTTTAAAAACGTAAAAAAACGTTTTTTTTTTTTCCCCT", "AAAAACGTAAAAAAACGT"),
                   0);
}

TEST(LevenshteinNonZero, LEVENSHTEIN_SUBSTRING) {
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("T", "A"), 1);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AA", "AT"), 1);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AA", "T"), 1);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AAAAAAAAACGGGGGGGCCCCCAAAAAAGGGGGGCCCCGGG", "T"), 1);

  EXPECT_DOUBLE_EQ(levenshtein_substring_score("AT", "TA"), 1);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("ATA", "TA"), 0);
  EXPECT_DOUBLE_EQ(levenshtein_substring_score("", "AA"), 2);
}
