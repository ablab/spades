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

void startSimpleGraphRecord(ostream &out, const string &name) {
	out << "digraph " << name << " {" << endl;
	out << "node" << "[";
	recordParameter(out, "fontname", "Courier");
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

double getColorParameter(int l, int r, double perc) {
	return l * perc + r * (1 - perc);
}

string getColor(int currentLength, int approximateLength) {
	currentLength %= approximateLength;
	int points[8][3] = {{0, 0, 1}, {0, 1, 1}, {1, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 1}, {0, 0, 1}};
	stringstream ss;
	int bound = approximateLength / 6;
	int num = currentLength / bound;
	double perc = (currentLength % bound) * 1. / bound;
	for(int i = 0; i < 3; i++) {
		ss << getColorParameter(points[num][i], points[num + 1][i], perc);
		if(i != 2)
			ss << ",";
	}
	return ss.str();
}

}
