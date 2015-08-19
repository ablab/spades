#include "err_helper_table.hpp"

#include <fstream>
#include <istream>

#include "logger/logger.hpp"

namespace hammer {
namespace errHelper {

namespace internal {

static const uint32_t helper_table_data[] = {
#include "err_helper_table.inc"
};

// numbers are cumulative sums of
// (2 * 4^^2) / 32,
// (2 * 4^^4) / 32,
// ...
const HelperTable helper_tables[] = {
  { 1, helper_table_data },
  { 2, helper_table_data + 1 },
  { 3, helper_table_data + 17 },
  { 4, helper_table_data + 273 },
  { 5, helper_table_data + 4369 }
};

}; // namespace internal

}; // namespace errHelper
}; // namespace hammer
