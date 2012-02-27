#include <stdinc.h>
#include <param.h>
#include <zigzag.h>
#include <extfunc.h>
#include <clean_graph.h>

int	noshave, filterinp;
int	qualinp;
char	inpfile[100], outfile[100], seqfile[100], qualfile[100];
char	htmlout, ***content, caption[2000];
FILE	*flog;

int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;


int main(int argc, char* argv[]) {

  if (argc != 7) {
    printf("Usage: simplifyGraph graphFile edgeFile intvFile readsFile vertexSize outBase\n");
    return 1;
  }

  char *graphFileName = argv[1];
  char *edgeFileName  = argv[2];
  char *intvFileName  = argv[3];
  char *readsFileName = argv[4];
  int vertexSize      = atoi(argv[5]);
  char *outputBase    = argv[6];

  NODES **vertices;
  EDGE  **edges;
  int numVertices, numEdges;
  READTABLE RT;
  /* Initialize A LOT! of constants */
  readpar();
  printf("Erosion length: %d\n", ErosionLength);

  read_graph_file(edgeFileName, graphFileName, &numVertices, &vertices, &numEdges, &edges);

  read_interval_file(intvFileName, numVertices, vertices);

  read_fasta_file(readsFileName, &RT);

  IGRAPH G;
  G.nodes = vertices;
  G.num_nodes = numVertices;
  

  CLEAN_PARAMS clean_params;

  /* parameters for cleaning the graph */
  clean_params.w = 5;
  clean_params.max_cov = 16;
  clean_params.C_B = BulgeLength;
  clean_params.L_b = BulgeCoverage;
  clean_params.C_W = WhirlLength;
  clean_params.L_c = ChimericCoverage;
  clean_params.L_e = 1;
  clean_params.do_zz = 1;              /* straighten zigzag paths */

  //	clean_params.T_c = ChimericTerm;
  //	clean_params.T_e = ErosionLength; 
  LABELVERT T_e_mem, T_c_mem;

  clean_params.T_c = &T_c_mem;
  if (ChimericTerm >= 0) {
    clean_params.T_c->depth = ChimericTerm;
    clean_params.T_c->dir = 0;
    clean_params.T_c->onepath = 1;     /* not used */
  } else {
    clean_params.T_c->depth = -ChimericTerm;
    clean_params.T_c->dir = 1;
    clean_params.T_c->onepath = 0;     /* not used */
  }
  clean_params.T_e = &T_e_mem;

  if (ErosionLength >= 0) {
    clean_params.T_e->depth = ErosionLength;
    clean_params.T_e->dir = 0;
    clean_params.T_e->onepath = 1;
  } else {
    clean_params.T_e->depth = -ErosionLength;
    clean_params.T_e->dir = 1;
    clean_params.T_e->onepath = 0;
  }

  clean_graph_once(&G, &RT, &clean_params);
  printf("read the graph file!!!\n");

  
  /* 
     set up the output file names.
  */
  char *edgeOutName = (char*) malloc(strlen(outputBase) + 6);
  strcpy(edgeOutName, outputBase);
  strcat(edgeOutName, ".edge");
  char *graphOutName = (char*) malloc(strlen(outputBase) + 7);
  strcpy(graphOutName, outputBase);
  strcat(graphOutName, ".graph");
  char *intvOutName = (char*) malloc(strlen(outputBase) + 6);
  strcpy(intvOutName, outputBase);
  strcat(intvOutName, ".intv");
  
  FILE *edgeOut, *graphOut, *intvOut;
  edgeOut = fopen(edgeOutName, "w");
  graphOut = fopen(graphOutName, "w");
  write_graph(vertices, numVertices, edgeOut, graphOut);
  fclose(edgeOut);
  fclose(graphOut);


  /* Output the intervals */
  write_interval_file(intvOutName, numVertices, vertices);


  char *gvzOutName = (char*) malloc(strlen(outputBase) + 5);
  strcpy(gvzOutName, outputBase);
  strcat(gvzOutName, ".dot");
  
  write_gvz_file(gvzOutName, numVertices, vertices, 0);
  
  return 0;
}
