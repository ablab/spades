#include <iostream>
#include <stdio.h>
#include "analyses.h"

int main (int argc,char *argv[])
{
	if  (argc <= 2) {
		std::cerr << "Не задан аргумент в командной строке.";
   		return 1;
 	}
	Analyses a(argv);
	return 0;
}
