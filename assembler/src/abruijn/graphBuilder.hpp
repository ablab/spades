#ifndef GRAPHBUILDER_H_
#define GRAPHBUILDER_H_

#include "graph.hpp"

class GraphBuilder
{
public:
	GraphBuilder() {}
	void build();
private:
    void selectGood();
};

#endif /* GRAPHBUILDER_H_ */
