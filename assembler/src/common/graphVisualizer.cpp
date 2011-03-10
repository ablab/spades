#include "graphVisualizer.hpp"

namespace gvis {

void startGraphRecord(ostream &out, const string &name) {
	out << "digraph " << name << " {" << endl;
}

void endGraphRecord(ostream &out) {
	out << "}" << endl;
}

string generateParameterString(const string &name, const string &value) {
	return name + "=" + "\"" + value + "\"";
}

}
