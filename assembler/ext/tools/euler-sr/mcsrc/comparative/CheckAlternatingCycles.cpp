/***************************************************************************
 * Title:          CheckAlternatingCycles.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <deque>
#include <unistd.h>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "InversionBins.h"
#include "mctypes.h"
#include "utils.h"

void GetColor(double weight, std::string &color);

void PrintGraphCharacter(IntVector &colors, 
			 std::ofstream &out, std::string &charTitle);




void BuildGraph(StringVector &species,
		InversionMatrix &invMatrix,
		IntMatrix &graph);

void PrintGraph(InversionMatrix &invMatrix,
		StringVector &species,
		IntMatrix &adjMat,
		IntVector &colors,
		FloatMatrix &scores,
		std::ofstream &out);

ssize_t ColorGraph(StringVector &species, 
	       IntMatrix &graph, 
	       IntVector &colors);

void AcyclicColor(IntMatrix &graph,
		  IntVector &colors);
		  

ssize_t blank, red, blue;
ssize_t same, opposite, unknown;

ssize_t GetToColor(ssize_t color) {
  if (color == red)
    return blue;
  if (color == blue)
    return red;
  return -1;
}

void PrintUsage() {
    std::cout << "usage: cacl invMatFile [-score scoresFile] [-graph graphVizFile]"
	      << "                    [-mat outputFileName] [-char charFileName] [-species speciesFile]  " << std::endl;
    std::cout << "                    [-i inversionNumber] [-concise] [-verbose verbosefile] " << std::endl;
    std::cout << "Look for odd cycles in a graph " << std::endl;
    std::cout << "       scoresFile   - a score matrix corresponding to the inversion matrix " << std::endl;
    std::cout << "       graphVizFile - a graphviz plot of the inversion graph will be put here " << std::endl;
    std::cout << "       charFileName - output the character as a colored traversal of th e grpah "<< std::endl;
    std::cout << "       inversionNumber - use this as the number for the inversion in the character file " << std::endl;
    std::cout << "       verbosefile  - print output here" << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    PrintUsage();
    exit(1);
  }
  std::string inversionFileName, 
    gvzFileName, 
    binOutFileName, 
    charFileName,
    scoresFileName,
    speciesFileName,
    verboseFileName;
  
  std::string inversionName;
  ssize_t fullCharOut     = 1;
  ssize_t i, j;
  ssize_t batchMode = 0;
  gvzFileName    = "";
  binOutFileName = "";
  charFileName   = "";
  scoresFileName  = "";
  speciesFileName = "";
  verboseFileName = "";
  inversionName = "";

  inversionFileName = argv[1];
  
  int argi = 2;
  argc--;
  ssize_t copt;
  std::ofstream vbout;
  while (argi <= argc) {
    if (strcmp(argv[argi], "-i") == 0) {
      inversionName = argv[++argi];
    }
    else if (strcmp(argv[argi], "-verbose") == 0) {
      verboseFileName = atoi(argv[++argi]);
      openck(verboseFileName, vbout, std::ios::out);
    }
    else if (strcmp(argv[argi], "-species") == 0) {
      speciesFileName = argv[++argi];
    }
    else if (strcmp(argv[argi], "-graph") == 0) {
      gvzFileName = argv[++argi];
    }
    else if (strcmp(argv[argi], "-mat") == 0) {
      binOutFileName = argv[++argi];
    }
    else if (strcmp(argv[argi], "-char") == 0) {
      charFileName = argv[++argi];
    }
    else if (strcmp(argv[argi], "-score") == 0) {
      scoresFileName = argv[++argi];
    }
    else if (strcmp(argv[argi], "-batch") == 0) {
      std::cout << "running in batch mode" << std::endl;
      batchMode = 1;
    }
    else if (strcmp(argv[argi], "-concise") == 0) {
      fullCharOut = 0;
    }
    else {
      std::cout << "did not understand option " << argv[argi] << std::endl;
      PrintUsage();
      return 1;
    }
    ++argi;
  }

  InversionList inversions;
  StringVector species;
  IntVector oddCycle;
  std::ofstream charOut;
  ReadInversionFile(inversionFileName, species, inversions);

  if (speciesFileName != "" ) {
    ReadSpeciesFile(speciesFileName, species);
  }
  if (charFileName != "") {
    openck(charFileName, charOut, std::ios::out);
    ssize_t s;
    if (fullCharOut ) {
      charOut << "species: ";
      for (s = 0; s < species.size() - 1; s++ ) {
	charOut << species[s] << " ";
      }
      charOut << species[s] << std::endl;
    }
  }

  // Read in auxillary information
  FloatMatrix scoreMatrix;
  if (scoresFileName != "") {
    std::ifstream scoresIn;
    openck(scoresFileName, scoresIn);
    ReadMatrix(scoresIn, scoreMatrix);
    // Normalize the score matrix;
    double max = -1;
    for (i = 0; i < scoreMatrix.size(); i++ ) {
      for (j = 0; j < scoreMatrix[i].size(); j++) {
	if (fabs(scoreMatrix[i][j]) > max) 
	  max = fabs(scoreMatrix[i][j]);
      }
    }
    ssize_t sign;
    for (i = 0; i < scoreMatrix.size(); i++ ) {
      for (j = 0; j < scoreMatrix[i].size(); j++) {
	sign = 1;
	if (scoreMatrix[i][j] < 0)
	  sign = -1;
	//	scoreMatrix[i][j] = (3.0/(4.0*max)) * scoreMatrix[i][j] + 0.25*sign;
	scoreMatrix[i][j] = scoreMatrix[i][j] / max;
      }
    }
    //    std::cout << "got matrix: " << std::endl;
    //    PrintMatrix(scoreMatrix, std::cout, scoreMatrix.size());
  }

  blank = 0; red = 1; blue = 2;
  same = 1; opposite = 0; unknown = 2;


  IntMatrix graph;
  ssize_t bin = 0;
  IntVector colors;
  bin = 0;
  ssize_t mat;
  for (mat = 0; mat < inversions.size(); ) {
    if (inversions[mat]->size() > 0) {
      BuildGraph(species, *inversions[mat], graph);
      ssize_t i, j;
      ssize_t foundOddCycle;
      if ((foundOddCycle = ColorGraph(species, graph, colors))) {
	// do something with the odd cycle.  Namely, remove this bin.
	oddCycle.push_back(mat);
      }

      if (foundOddCycle) {
	std::cout << "odd cycle found for inversion " << inversions[mat]->title << std::endl;
      }
      if (!foundOddCycle and charFileName != "") {
	//	AcyclicColor(graph, colors);
	if (inversions[mat]->title != "") 
	  PrintGraphCharacter(colors, charOut, inversions[mat]->title);
	else {
	  if (inversionName == "") { 
	    std::stringstream namestrm;
	    namestrm << mat;
	    inversionName = namestrm.str();
	  }
	  PrintGraphCharacter(colors, charOut, inversionName);
	}
      }
      if (gvzFileName != "") {
	std::cout << inversions[mat]->title << std::endl;
	std::cout << "writing to " << gvzFileName << " must enter a char to continue ";
	std::ofstream gvzOut;
	openck(gvzFileName, gvzOut, std::ios::out);
	std::cout << scoreMatrix.size() << std::endl;
	PrintGraph(*inversions[mat],
		   species, graph, colors,
		   scoreMatrix, gvzOut);
	gvzOut.close();
        std::cout << "done writing " << std::endl;
	char c;
	if (! batchMode) 
	  c = std::cin.get();
	if (!batchMode and c == 'p' and mat != 0) {
	  std::cin.get(); // get the eol
	  --mat;
	}
	else {
	  ++mat;
	}
      }
      else {
	++mat;
      }
    }
  }
  if (binOutFileName != "") {
    std::ofstream binOut;
    // Remove bins that corresponded to odd cycl
    std::sort(oddCycle.begin(), oddCycle.end());
    mat = inversions.size()-1; 
    ssize_t cycleIdx = oddCycle.size() - 1;
    while ( mat >= 0 and cycleIdx >= 0) {
      if (mat == oddCycle[cycleIdx]) {
	inversions.erase(inversions.begin() + mat);
	cycleIdx--;
      }
      else
	--mat;
    }
    openck(binOutFileName, binOut, std::ios::out);
    PrintBins(species, inversions, binOut, 0);
    std::cout.flush();
    binOut.close();
  }
  return 0;
}

ssize_t ColorGraph(StringVector &species, IntMatrix &graph, IntVector &colors) {
  
  ssize_t V;
  V = graph.size();
  colors.resize(graph.size());
  ssize_t i, j;
  std::set<ssize_t> vertices;
  std::deque<ssize_t> toVisit;
  for (i = 0; i < V; i++) {
    colors[i] = blank;
    vertices.insert(i);
  }
  IntMatrix connectedComponents;
  ssize_t ncc = 0;
  ssize_t edgeFound;
  ssize_t u, v;
  // Vertices contains a list of unvisited nodes.
  while (vertices.size() > 0) {
    edgeFound = 0;
    while (!edgeFound and vertices.size() > 0) {
      u = *(vertices.begin());
      v = 0;
      while (!edgeFound and v < V) {
	if (graph[u][v] == opposite) {
	  edgeFound = 1;
	}
	else {
	  v++;
	}
      }
      if (!edgeFound ) {
	// vertex u is not connected to anything
	// color it so that it is not looked at agin
	// and remove it
	vertices.erase(u);
      }
    }
    if (!edgeFound) 
      // No edges left in the graph
      break;
    ssize_t fromColor, toColor;
    colors[u] = red;
    toColor = GetToColor(colors[u]);
    toVisit.push_back(u);
    ncc++;
    connectedComponents.resize(ncc);

    while (toVisit.size() > 0) {
      u = toVisit.front(); toVisit.pop_front();
      connectedComponents[ncc-1].push_back(u);
      fromColor = colors[u];
      toColor   = GetToColor(fromColor);
      vertices.erase(u);  
      // Check for the edges from this vertex
      for (v = 0; v < V; v++) {
	// Don't consider self edges 
	// they shouldn't exist anyway, so maybe I should keep the check?
	if (v == u or graph[u][v] == same or graph[u][v] == unknown) continue;
	assert(graph[u][v] == opposite);
	if (colors[v] != blank) {
	  if (colors[v] == fromColor) {
	    // the to color should always be different from from color
	    return 1;
	  }
	}
	else {
	  colors[v] = toColor;
	  toVisit.push_back(v);
	}
      }
    }
  }
  return 0; 
}    
    
  
void BuildGraph(StringVector &species,
		InversionMatrix &invMatrix,
		IntMatrix &adjMat) {

  CreateMatrix(adjMat, species.size(), species.size());
  ssize_t r, c;
  // create an n x n adjacency matrix
  // from m x n input data, where the 2*(n-m) unknown
  // values remain as unknown
  for (r = 0; r < species.size(); r++ )
    for (c = 0; c < species.size(); c++ )
      adjMat[r][c] = unknown;
  
  ssize_t inv, spec, spec2;
  for (inv = 0; inv < invMatrix.size(); inv++) {
    spec = FindSpecies(species, invMatrix.loci[inv]->species);
    if (spec != -1) {
      for (spec2 = 0; spec2 < invMatrix.loci[inv]->size(); ++spec2) {
	if ( invMatrix.loci[inv]->orient[spec2] == 0 ) {
	  adjMat[spec][spec2] = opposite;
	  adjMat[spec2][spec] = opposite;
	}
	else if (invMatrix.loci[inv]->orient[spec2] == 1) {
	  adjMat[spec2][spec] = same;
	  adjMat[spec][spec2] = same;
	}
	else {
	  adjMat[spec][spec2] = unknown;
	  adjMat[spec2][spec] = unknown;
	}
      }
    }
  }
}

void PrintGraphCharacter(IntVector &colors, 
			 std::ofstream &out, std::string &binName) {

  ssize_t i;

  // Find the maximum color.
  ssize_t nr, nb;
  nr = nb = 0;
  for (i = 0; i < colors.size(); i++ ){
    if (colors[i] == red)
      nr++;
    else if (colors[i] == blue)
      nb++;
  }
  ssize_t redNumber, blueNumber;
  if (nr > nb) {
    redNumber = 1;
    blueNumber = 0;
  }
  else {
    redNumber = 0;
    blueNumber = 1;
  }

  out << binName << "\t";
  for (i = 0; i < colors.size(); i++ ){
    if (colors[i] == blank) {
      out << " " << 2;
    }
    else if (colors[i] == red) {
      out << " " << redNumber;
    }
    else if (colors[i] == blue) {
      out << " " << blueNumber;
    }
  }
  out << std::endl;
}


void PrintGraph(InversionMatrix &invMatrix,
		StringVector &species,
		IntMatrix &adjMat,
		IntVector &colors,
		FloatMatrix &scores,
		std::ofstream &out) {
  out << "digraph G { " << std::endl;
  out << "size=\"8,10\";" << std::endl;
  ssize_t u,v;
  // Print the first group
  std::vector<ssize_t> redIndex, blueIndex;
  out << "subgraph clusterBLUE { " << std::endl;
  out << "color = white;" << std::endl;

  for (u = 0; u < species.size(); u++) {
    if (colors[u] == blue) 
      blueIndex.push_back(u);
    if (colors[u] == red)
      redIndex.push_back(u);
  }
  ssize_t i;
  for (i = 0; i < blueIndex.size(); i++ ) {
    u = blueIndex[i];
    out << u << " [ label=\""<< species[u] << "\"]" << std::endl;
    out << u << " [ fillcolor=\"blue\"]" << std::endl;
    out << u << " [style=\"filled\"]" << std::endl;
    out << u << " [width=\"1.5\"] " << std::endl;
    out << u << " [pos=\""<<0 << "," << u*2 << "!\"]" << std::endl;
    if (i < blueIndex.size()-1) {
      out << blueIndex[i] << " -> " << blueIndex[i+1] 
	  << " [style=\"invis\"][weight=\"1000\"]; " << std::endl;
    }
  }
  out << "} " << std::endl;

  out << "subgraph clusterRED { " << std::endl;
  out << "color = white;" << std::endl;
  for (i = 0; i < redIndex.size(); i++ ) {
    u = redIndex[i];
    out << u << " [ label=\""<< species[u] << "\"]" << std::endl;
    out << u << " [ fillcolor=\"red\"]" << std::endl;
    out << u << " [style=\"filled\"]" << std::endl;
    out << u << " [pos=\""<<2 << "," << u*2 << "!\"]" << std::endl;
    out << u << " [width=\"1.5\"] " << std::endl;
    out << u << " [fixedsize=\"true\"] " << std::endl;
    if (i < redIndex.size()-1) {
      out << redIndex[i] << " -> " << redIndex[i+1] 
	  << " [style=\"invis\"][weight=\"1000\"]; " << std::endl;
    }
  }
  out << "} " << std::endl;
  for (u = 0; u < species.size(); u++) {
    if (colors[u] == blank) {
      out << u << " [ label=\""<< species[u] << "\"]" << std::endl;
      out << u << " [ fillcolor=\"gray\"]" << std::endl;
      out << u << " [style=\"filled\"]" << std::endl;
      out << u << " [width=\"1.5\"] " << std::endl;
      out << u << " [fixedsize=\"true\"] " << std::endl;
    }
  }
  ssize_t useEdgeColoring = 0;
  if (scores.size() > 0) {
    std::cout << "using edge coloring " << std::endl;
    useEdgeColoring = 1;
  }
  ssize_t from, to;
  std::string colorStr;

  for (u = 0; u < species.size() -1 ; u++) {
    for (v = u+1; v < species.size(); v++ ) { 
      from = FindSpecies(invMatrix, species[u]);
      to   = FindSpecies(invMatrix, species[v]);
      if (adjMat[u][v] != unknown) {
	if ((useEdgeColoring and from != -1 and to != -1)
	    or adjMat[u][v] == opposite or adjMat[v][u] == opposite) {
	  out << u << " -> " << v;
	}
	if (useEdgeColoring and from != -1 and to != -1) {
	  GetColor(scores[from][to], colorStr);
	  out << " [ color = \"" << colorStr << "\" ] ";
	}
	if ((useEdgeColoring and from != -1 and to != -1)
	    or adjMat[u][v] == opposite) 
	  out << " [style=\"setlinewidth(2)\"] [samehead][sametail] [arrowhead = \"none\"];" 
	      << std::endl;
      }
    }
  }
  out << " }" << std::endl;
}


void GetColor(double weight, std::string &color) {
  ssize_t blueMat[48] = {   0,0,255, 0,0,238, 0,0,221, 0,0,204, 
			0,0,187, 0,0,170, 0,0,153, 0,0,136,
			0,0,119, 0,0,102, 0,0,85 , 0,0,68 ,
			0,0,51 , 0,0,34 , 0,0,17 , 0,0,0 };

  ssize_t redMat[48] = {255, 0,0, 238, 0,0, 221, 0,0, 204, 0,0, 
		    187, 0,0, 170, 0,0, 153, 0,0, 136, 0,0,    
		    119, 0,0, 102, 0,0, 85 , 0,0, 68 , 0,0,    
		    51 , 0,0, 34 , 0,0, 17 , 0,0, 0  , 0,0 };

  if (weight > 0) {
    std::stringstream colorStream;
    // Create a blue edge
    // Blue coloring
    ssize_t rowIndex = (ssize_t) std::floor(((1 - weight) * 16));
    colorStream << "0.666," << weight  << ",1";
    /*std::ios::hex << (1.0*blueMat[rowIndex*3])
	<< std::ios::hex << (1.0*blueMat[rowIndex*3 + 1])
    colorStream <<  blueMat[rowIndex*3 + 2];
    colorStream.setf(std::ios_base::dec, std::ios_base::basefield);
    */
    color = colorStream.str();
  }
  else {
    std::stringstream colorStream;
    weight = -weight;
    colorStream << "0," << weight << ",1";
    /*    colorStream.setf(std::ios_base::hex, std::ios_base::basefield);
    colorStream << (redMat[rowIndex*3]) << "0000";
		<< std::ios::hex << (1.0*redMat[rowIndex*3 + 1])
		<< std::ios::hex << (1.0*redMat[rowIndex*3 + 2]);
    */
    color = colorStream.str();
  }
}

void AcyclicColor(IntMatrix &graph,
		  IntVector &colors) {
  // It is possible that several vertices are not 
  // along cycles, but are connected to vertices 
  // that are colored, but in the same orientation.

  // Check that here.
  ssize_t u, v;
  for (u = 0; u < colors.size(); u++ ) {
    if (colors[u] == blank) {
      for (v = 0; v < colors.size(); v++ ) {
	if (u == v or colors[v] == blank) continue;
	if (graph[u][v] == same) {
	  colors[u] = colors[v];
	  break;
	}
      }
    }
    if (colors[u] == blank) {
      for (v = colors.size()-1; v >= 0; --v ) {
	if (u == v or colors[v] == blank) continue;
	if (graph[u][v] == same) {
	  colors[u] = colors[v];
	  break;
	}
      }
    }
  }
}
