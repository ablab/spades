#include "err_helper_table.hpp"

#include <fstream>
#include <istream>

#include "logger/logger.hpp"

namespace hammer {
namespace errHelper {

void initHelperTables(const std::string& filename) {
  std::ifstream stream(filename);
  for (size_t i = 1; i <= internal::MAX_K; ++i)
    internal::helper_tables.push_back(internal::HelperTable(i, stream));
}

namespace internal {

std::vector<HelperTable> helper_tables;

HelperTable::HelperTable(unsigned k, std::istream& stream) {
  size_t sz = 1 << (4 * k); // # of hints - 4 ^^ 2k 
  size_t bytes = sz / 4; // each char contains 4 hints
  storage_.resize(bytes);

  VERIFY(stream.good());
  stream.read(&storage_[0], bytes);
}

}; // namespace internal

}; // namespace errHelper
}; // namespace hammer
