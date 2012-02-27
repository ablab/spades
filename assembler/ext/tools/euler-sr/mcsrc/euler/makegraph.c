/***************************************************************************
 * Title:          makegraph.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

void makegraph(IGRAPH *G, READTABLE *RT, ALIGN **eq_class);
void destroygraph(IGRAPH *G);

void makegraph(IGRAPH *G, READTABLE *RT, ALIGN **eq_class)
{
  make_Agraph(RT, eq_class, 1, 1, 0,
	      G);
}

/***********************************************************************
 * make_Agraph(RT, eq_class, split_badnodes, is_sym, use_int,
 *             G)
 *
 * input:
 *   RT->num_seq, src_seq, len_seq: table of reads
 *   eq_class: alignments
 *   split_badnodes:
 *        1 = split bad nodes
 *        0 = don't
 *   is_sym: NOT YET IMPLEMENTED
 *        1 = the reads, nodes, and edges come in complementary pairs
 *            (CURRENTLY IT WORKS THIS WAY)
 *
 *        0 = they don't
 *            (NOT YET IMPLEMENTED)
 *
 *   use_int: NOT YET IMPLEMENTED
 *        1 = produce interval graph
 *        0 = produce nucleotide graph (CURRENTLY IT WORKS THIS WAY)
 *
 * output:
 *   G: graph
 *   RT->readlist: how reads map into the graph
 ***********************************************************************/

void make_Agraph(READTABLE *RT,
		 ALIGN **eq_class,
		 int split_badnodes,
		 int is_sym,   /* reserved: 1=symmetric graph, 0=not */
		 int use_int,  /* reserved: 1=interval graph, 0=nuc graph */

		 IGRAPH *G
		 )
{
  int	i, j, k, l, m, n, m0, n0;
  int	num_seq, num_edge, *len_seq;
  char	**src_seq;
  READLIST	**readlist;
  int	num_vertex, nbad;
  NODES	**vertex, **badnode;

  num_seq = RT->num_seq;
  src_seq = RT->src_seq;
  len_seq = RT->len_seq;

  /*	Initialize the nodes: each position in each read is assigned
	as a new node. An array of "readlist" is set up for each read	*/

  readlist = (READLIST **) ckalloc(2 * num_seq * sizeof(READLIST *));
  for (i = 0; i < 2 * num_seq; i ++)	{
    readlist[i] = (READLIST *) ckalloc(len_seq[i] * sizeof(READLIST));
  }
  printf("intitialize nodes...\n");
  initialize(readlist, len_seq, num_seq);
  printf("done.\n");
  num_vertex = countnode(readlist, len_seq, 2 * num_seq);
  printf("# of nodes before merge: %d\n", num_vertex);
  //  cleannode(readlist, len_seq, num_seq);

  /*	Glue together two nodes if their corresponding positions are defined
	as equivalent in a pairwise alignment.		*/

  printf("Merge...\n");
  merge_gap(src_seq, num_seq, len_seq, eq_class, readlist, gap_k);
  printf("done.\n");

  /* Count the new number of unique vertices */
  num_vertex = countnode(readlist, len_seq, 2 * num_seq);
  printf("# of nodes after merge: %d\n", num_vertex);
  cleannode(readlist, len_seq, num_seq);


  /* Mark and fix bad supernodes */
  if (split_badnodes) {

    /*	Check if one super node is bad or not	*/
    badnode = (NODES **) ckalloc(num_vertex * sizeof(NODES *));
    nbad = 0;
    for (i = 0; i < 2 * num_seq; i ++)	{
      for (j = 0; j < len_seq[i]; j ++)	{
	if(!readlist[i][j].node->visit)	{
	  readlist[i][j].node->num_path = countthickness(readlist[i][j].node);
	  if(readlist[i][j].node->num_path > 1)	{
	    badnode[nbad ++] = readlist[i][j].node;
	  }
	  readlist[i][j].node->visit = 1;
	}
      }
    }
    cleannode(readlist, len_seq, num_seq);

  /*	Create separate supernodes for each positions in the supernodes	*/

    printf("separating...\n");
    for (k = 5; k >= 1; k --)	{
      m = n = nbad + 1;
      do	{
	m0 = m;
	n0 = n;
	m = n = 0;
	for (i = 0; i < nbad; i ++)	{
	  if(badnode[i] && badnode[i]->num_path > 1)	{

	    /*	Algorithm 0 for spliting bad supernodes		*/
	    /*
	      separatenode(badnode[i], readlist, len_seq, num_seq);
	    */
	    /*	Algorithm 1 for spliting bad supernodes		*/
	    separatenode0(badnode[i], readlist, len_seq, num_seq, k);
	    n ++;
	    m += badnode[i]->npos;
	  }
	}
	l = 0;
	for (i = 0; i < nbad; i ++)	{
	  if(badnode[i] && badnode[i]->npos == 0)	{
	    badnode[i] = free_nodes(badnode[i]);
	  }
	  if(badnode[i] && badnode[i]->num_path > 1)	{
	    badnode[i]->num_path = countthickness(badnode[i]);
	    if(badnode[i]->num_path > 1)	{
	      l ++;
	    }
	  }
	}
      } while(n < n0 || m < m0);
      printf("%d bad supernodes (%d read positions) separated with threshold %d.\n", n, m, k);
      printf("%d bad supernodes remain.\n", l);
    }
    for (i = 0; i < nbad; i ++)	{
	  if(badnode[i] && badnode[i]->num_path > 1)	{
	      separatenode(badnode[i], readlist, len_seq, num_seq);
	  }
	  if(badnode[i] && badnode[i]->npos == 0)	{
	    badnode[i] = free_nodes(badnode[i]);
	  }
    }
    free((void **) badnode);
    printf("done.\n");

    /*	Count # of nodes and allocate the memory of the pointer to these nodes	*/
    
    num_vertex = countnode(readlist, len_seq, 2 * num_seq);
    printf("# of nodes after separating: %d\n", num_vertex);
    cleannode(readlist, len_seq, num_seq);
  }

  vertex = (NODES **) ckalloc(num_vertex * sizeof(NODES *));

  /*	Build nucleotide graph	*/

  printf("Creating edges...\n");
    fflush(stdout);

    /* vertex is a list of unique nodes (should be of length 'num_vertex' above).
       the vertices are connected if there is a sequence that contains two adjacent nucleotides
       corresponding to two vertices in the vertex list.
    */
  num_edge = link_edge(num_seq, len_seq, readlist, vertex);
  printf("done.\n");
    fflush(stdout);

  /*	Check the number of edges from a vertex to itself (self-loops)	*/

  n = 0;
  for (i = 0; i < num_vertex; i ++)	{
    for (j = 0; j < vertex[i]->num_nextedge; j ++)	{
      if(vertex[i]->nextedge[j]->end == vertex[i])	{
	n ++;
      }
    }
  }
  printf("# of self-loops: %d.\n", n);
    fflush(stdout);

    /* 
       Store everything in an easily accessible list structure to be used by glenn's code.
    */
  RT->readlist = readlist;
  G->nodes = vertex;
  G->num_nodes = num_vertex;
}

void destroygraph(IGRAPH *G)
{
  free_graph(G->nodes, G->num_nodes);
  free((void **) G->nodes);
}
