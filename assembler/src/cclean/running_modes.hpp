#ifndef RUNNING_MODES_H_
#define RUNNING_MODES_H_

#include <ostream>
#include <string>
#include "additional.cpp"

class ireadstream;
namespace cclean {
class AdapterIndex;
}

using additional::WORK_MODE_TYPE;

void exactAndAlign(std::ostream &output, std::ostream &bed, ireadstream *input,
                   const std::string &db, const cclean::AdapterIndex &index,
                   const WORK_MODE_TYPE &brute);

#endif /* RUNNING_MODES_H_ */
