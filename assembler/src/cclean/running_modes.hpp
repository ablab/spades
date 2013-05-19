#ifndef RUNNING_MODES_H_
#define RUNNING_MODES_H_

#include <ostream>

class ireadstream;
namespace cclean {
class AdapterIndex;
}

void exactAndAlign(std::ostream& output, std::ostream& bed, ireadstream * input, const cclean::AdapterIndex &index);

#endif /* RUNNING_MODES_H_ */
