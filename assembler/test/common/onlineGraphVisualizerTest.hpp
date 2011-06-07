#include "graphVisualizer.hpp"
#include "iostream"
#include "string"
#include "cute.h"

using namespace gvis;

void testOnline() {
	stringstream ss;
	PairedGraphPrinter<int> g1("", ss);
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph  {\nnode[fontname=<Courier>,shape=<plaintext>]\n", tmp);
}

void testOnlineEmpty() {
	stringstream ss;
	PairedGraphPrinter<int> g("", ss);
	g.open();
	g.close();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph  {\nnode[fontname=<Courier>,shape=<plaintext>]\n}\n", tmp);
}

void testOnlineEmptyWithName() {
	stringstream ss;
	PairedGraphPrinter<int> g("myName", ss);
	g.open();
	g.close();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph myName {\nnode[fontname=<Courier>,shape=<plaintext>]\n}\n", tmp);
}

void testOnlineSingleVertex() {
	stringstream ss;
	DotGraphPrinter<int> g ("myName", ss);
	g.open();
	g.AddVertex(0, "oppa");
	g.close();
	string tmp;
	//	cout <<
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph myName {\nnode[fontname=<Courier>]\nvertex_0[label=<oppa>,style=<filled>,color=<black>,fillcolor=<white>]\n}\n", tmp);
}

void testOnlineSingleVertexPair() {
	stringstream ss;
	PairedGraphPrinter<int> g("myName", ss);
	g.open();
	g.AddVertex(1, "AAAAAAAAA", 2, "TTTTTTTTT");
	g.close();
	string tmp = ss.str();
	ASSERT_EQUAL("digraph myName {\nnode[fontname=<Courier>,shape=<plaintext>]\nvertex_1_2[label=<<TABLE>\n<TR><TD BORDER = \"0\" PORT = \"port_1_in\"></TD><TD BORDER = \"0\" PORT = \"port_\">AAAAAAAAA</TD><TD BORDER = \"0\" PORT = \"port_1_out\"></TD></TR>\n<TR><TD BORDER = \"0\" PORT = \"port_2_out\"></TD><TD BORDER = \"0\" PORT = \"port_\">TTTTTTTTT</TD><TD BORDER = \"0\" PORT = \"port_2_in\"></TD></TR>\n</TABLE>>,style=<filled>,color=<black>,fillcolor=<white>]\n}\n", tmp);
}

void testOnlineSingleComplexEdge() {
	stringstream ss;
	PairedGraphPrinter<int> g("myName", ss);
	g.open();
	g.AddVertex(1, "AAAAAAAAA", 2, "TTTTTTTTT");
	g.AddVertex(3, "CCCCCCCCC", 4, "GGGGGGGGG");
	g.AddEdge(make_pair(2,1), make_pair(3,4), "oppa");
	g.close();
	string tmp = ss.str();
	ASSERT_EQUAL("digraph myName {\nnode[fontname=<Courier>,shape=<plaintext>]\nvertex_1_2[label=<<TABLE>\n<TR><TD BORDER = \"0\" PORT = \"port_1_in\"></TD><TD BORDER = \"0\" PORT = \"port_\">AAAAAAAAA</TD><TD BORDER = \"0\" PORT = \"port_1_out\"></TD></TR>\n<TR><TD BORDER = \"0\" PORT = \"port_2_out\"></TD><TD BORDER = \"0\" PORT = \"port_\">TTTTTTTTT</TD><TD BORDER = \"0\" PORT = \"port_2_in\"></TD></TR>\n</TABLE>>,style=<filled>,color=<black>,fillcolor=<white>]\nvertex_3_4[label=<<TABLE>\n<TR><TD BORDER = \"0\" PORT = \"port_3_in\"></TD><TD BORDER = \"0\" PORT = \"port_\">CCCCCCCCC</TD><TD BORDER = \"0\" PORT = \"port_3_out\"></TD></TR>\n<TR><TD BORDER = \"0\" PORT = \"port_4_out\"></TD><TD BORDER = \"0\" PORT = \"port_\">GGGGGGGGG</TD><TD BORDER = \"0\" PORT = \"port_4_in\"></TD></TR>\n</TABLE>>,style=<filled>,color=<black>,fillcolor=<white>]\nvertex_1_2:port_2_out->vertex_3_4:port_3_in[label=<oppa>,color=<black>]\n}\n", tmp);
}

void testOnlineSingleEdge() {
	stringstream ss;
	DotGraphPrinter<int> g("myName", ss);
	g.open();
	g.AddEdge(0, 1, "oppa");
	g.close();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph myName {\nnode[fontname=<Courier>]\nvertex_0->vertex_1[label=<oppa>,color=<black>]\n}\n", tmp);
}


//void testOnlineFile() {
//	PairedGraphPrinter<int> *g = new PairedGraphPrinter<int> ("", "test/data/oppa.txt");
//	g->addVertex(0, "oppa");
//	g->addEdge(0, 1, "oppa");
//	g->output();
//
//}

cute::suite onlineGraphVisualizerSuite() {
	cute::suite s;
	s.push_back(CUTE(testOnline));
	s.push_back(CUTE(testOnlineEmpty));
	s.push_back(CUTE(testOnlineEmptyWithName));
	s.push_back(CUTE(testOnlineSingleVertex));
	s.push_back(CUTE(testOnlineSingleEdge));
//	s.push_back(CUTE(testOnlineFile));
	s.push_back(CUTE(testOnlineSingleVertexPair));
	s.push_back(CUTE(testOnlineSingleComplexEdge));
	return s;
}
