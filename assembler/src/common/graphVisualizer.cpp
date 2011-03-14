#include "graphVisualizer.hpp"

namespace gvis {

void startGraphRecord(ostream &out, const string &name) {
	out << "digraph " << name << " {" << endl;
	out << "node" << "[";
	recordParameter(out, "fontname", "Courier");
	out << ",";
	recordParameter(out, "shape", "plaintext");
	out << "]" << endl;
}

void endGraphRecord(ostream &out) {
	out << "}" << endl;
}

void recordParameter(ostream &out, const string &name, const string &value) {
	out << name << "=" << "<" << value << ">";
}

string constructCell(const string &label, int border, const string &port) {
	stringstream ss;
	ss << "<TD BORDER = \"" << border << "\" PORT = \"port_" << port << "\">"
			<< label << "</TD>";
	return ss.str();
}

}
