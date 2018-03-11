#include <gtest/gtest.h>

#include "graph.hpp"
#include "demo.hpp"
#include "fees.hpp"

double levenshtein_substring_score(const std::string &s, const std::string &query) {
  auto fees = hmm::levenshtein_fees(query);
  auto graph = Graph({s});
  auto result = find_best_path(fees, graph.all());
  auto paths = result.top_k(100);  
  return paths[0].second;
}

TEST(LevenshteinZero, LEVENSHTEIN_SUBSTRING) {
  EXPECT_FLOAT_EQ(levenshtein_substring_score("A", "A"), 0);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("AAAAA", "A"), 0);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("AAAAACGTAAAAAAACGT", "AAA"), 0);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("AAAAACGTAAAAAAACGT", ""), 0);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("", ""), 0);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("AAAAACGTAAAAAAACGT", "CGT"), 0);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("AAAAACGTAAAAAAACGT", "AAAAACGTAAAAAAACGT"), 0);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("TTTTTTTTTTTTTAAAAACGTAAAAAAACGTTTTTTTTTTTTCCCCT", "AAAAACGTAAAAAAACGT"),
                  0);
}

TEST(LevenshteinNonZero, LEVENSHTEIN_SUBSTRING) {
  EXPECT_FLOAT_EQ(levenshtein_substring_score("T", "A"), 1);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("AA", "AT"), 1);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("AA", "T"), 1);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("AAAAAAAAACGGGGGGGCCCCCAAAAAAGGGGGGCCCCGGG", "T"), 1);

  EXPECT_FLOAT_EQ(levenshtein_substring_score("AT", "TA"), 1);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("ATA", "TA"), 0);
  EXPECT_FLOAT_EQ(levenshtein_substring_score("", "AA"), 2);
}
