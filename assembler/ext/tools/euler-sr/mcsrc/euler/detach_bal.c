/***************************************************************************
 * Title:          detach_bal.c
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
#include <assert.h>

void CollectPaths(EDGE *inEdge, EDGE *outEdge, NODES *vertex, PATH *paths, int numPaths, 
									int**collectedPaths, int *numCollectedPaths) {

	int vp, vpe;
	int pathIndex, pathEdgeIndex;
	*numCollectedPaths = 0;
	for (vp = 0; vp < vertex->num_path; vp++) {
		pathIndex = vertex->path_index[vp];
		pathEdgeIndex  = vertex->path_pos[vp];
		if (pathEdgeIndex > 0 and
				paths[pathIndex].len_path > 1) {
			if (paths[pathIndex].edge[pathEdgeIndex-1] == inEdge and
					paths[pathIndex].edge[pathEdgeIndex] == outEdge) {
				(*numCollectedPaths)++;
			}
		}
	}
	*collectedPaths = (int*) ckalloc(sizeof(int) * (*numCollectedPaths));
	int pathNumber = 0;
	for (vp = 0; vp < vertex->num_path; vp++) {
		pathIndex = vertex->path_index[vp];
		pathEdgeIndex  = vertex->path_pos[vp];
		if (pathEdgeIndex > 0 and
				paths[pathIndex].len_path > 1) {
			if (paths[pathIndex].edge[pathEdgeIndex-1] == inEdge and
					paths[pathIndex].edge[pathEdgeIndex] == outEdge) {
				(*collectedPaths)[pathNumber] = pathIndex;
				++pathNumber;
			}
		}
	}
}


void PrintIndexedPaths(PATH *paths, int *pathIndices, int numIndices){ 
	int index, pos;
	for (index = 0; index < numIndices; index++) {
		printf("path: %d ", pathIndices[index]);
		for (pos = 0; pos < paths[pathIndices[index]].len_path; pos++) {
			printf(" %d ", paths[pathIndices[index]].edge[pos]->index);
		}
		printf("\n");
	}
}

void PrintPaths(PATH *paths, int numPaths) {
	int path, pos;
	for (path = 0; path < numPaths; path++ ){ 
		printf("path: %d ", path);
		for (pos = 0; pos < paths[path].len_path; pos++) {
			printf(" %d ", paths[path].edge[pos]->index);
		}
		printf("\n");
	}
}

EDGE *detach_bal(EDGE *edge1, EDGE *edge2, PATH *path, int num_path, int *f_edge, int num_seq);

EDGE *detach_bal(EDGE *edge1, EDGE *edge2, PATH *path, int num_path, int *f_edge, int num_seq)
{
  /* Input:
     edge1 and edge2 are connected through a vertex:
     edge1 --> vertex --> edge2

     f_edge boolean array indicating whether edge1 or edge2 are freed. 
  */
  
  /*  There exists --edge1-->vertex--edge2--> in the graph, and it should be ok to do equivalent
      transformation on them.  All paths that go through edge1 to edge2 should be re-routed 
      so that they don't use 'vertex'.
  */
	int	i, j, k, l, m, n, q, n1, n2, c, min_c, min_l, min_i, num_r, l1, l2;
	double	r, min_r, min_min_r;
	int	num_endpath, num_startpath, num_midforpath, num_midaftpath;
	int	num_startmatch, num_endmatch;
	int	*num_match1, *num_match2, *num_cont;
	int	*edgematch1, *edgematch2;
	PATH	*endpath, *startpath;
	PATH	*midforpath, *midaftpath;
	NODES	*begin, *end, *vertex;
	int	*beginlist, *endlist, nbegin, nend;
	EDGE	*edge, *newedge;

	if(edge1 -> end != edge2 -> begin)	{
		printf("edge1 %d (%d) %d edge2 %d (%d) %d\n", edge1->index, edge1, edge1 -> end, 
					 edge2->index, edge2, edge2 -> begin);
		printf("Edge not connected\n");
		return(NULL);
	}

	f_edge[0] = f_edge[1] = 0;

	vertex = edge1 -> end;

	int *pathIndices, numPaths;
	/*
		CollectPaths(edge1, edge2, vertex, path, num_path, &pathIndices, &numPaths);
		printf("the paths through %d -> %d -> %d are: \n", edge1->index, vertex->index, edge2->index);
		PrintIndexedPaths(path, pathIndices, numPaths);
		
		free((void*) pathIndices);
	*/

/*	num_startmatch--
	num_endmatch  --number of branches consistent with edge2	*/

	
	num_startmatch = countstartmatch(edge1, vertex, path, num_path);
	num_endmatch   = countendmatch(edge2, vertex, path, num_path);
	/*
		printf("  detaching edge %d, edge %d paths starting e1: %d (len %d) paths starting e2: %d (len %d)\n", 
				 edge1->index, edge1->length, edge2->index, edge2->length, num_startmatch, num_endmatch);


				 printf("edge: %d (bal: %d) %d -> %d num_startmatch: %d  edge %d (bal: %d) %d -> %d num_endmatch: %d\n",
	       edge1, edge1->bal_edge, edge1->begin->index, edge1->end->index, num_startmatch,
	       edge2, edge2->bal_edge, edge2->begin->index, edge2->end->index, num_endmatch);
	*/
	if(num_startmatch == 0 || num_endmatch == 0)	{
		printf("Path not found %d(%d,%d-%d) %d %d(%d,%d-%d) %d \n", edge1, edge1 -> length,
			edge1 -> begin, edge1 -> end,
			num_startmatch, edge2, edge2 -> length, edge2 -> begin, edge2 -> end, num_endmatch);
		printf("vertex %d %d %d num_path %d\n", vertex, vertex -> num_lastedge, vertex -> num_nextedge,
			vertex -> num_path);
		/*		assert(0);*/
		return NULL;
	}

	endpath   = (PATH *) ckalloc(vertex -> num_path * sizeof(PATH));
	startpath = (PATH *) ckalloc(vertex -> num_path * sizeof(PATH));
	/* Store the paths that end in edge1. */
	num_endpath = collectendpaths(vertex, edge1, endpath, path, num_path); /* P->x */
	/*
		printf("endpaths edge: %d (len %d): \n", edge1->index, edge1->length);
		PrintPaths(endpath, num_endpath);
	*/

	/* Collect paths that start in edge2. */
	num_startpath = collectstartpaths(vertex, edge2, startpath, path, num_path); /* Py-> */

	midforpath = (PATH *) ckalloc(vertex -> num_path * sizeof(PATH));
	midaftpath = (PATH *) ckalloc(vertex -> num_path * sizeof(PATH));
	num_cont = (int *) ckalloc((num_startpath + num_endpath) * sizeof(int));
	for(j = 0; j < num_endpath; j ++) num_cont[j] = 0;


/*	num_match1--number of branches consistend with P->x	*/

	num_match1 = (int *) ckalloc(num_endpath * sizeof(int));

	/* edgematch1 - a boolean array that is 1 if the endpath 'i' is consistent
	                with just 1 out edge from 'vertex', and 0 otherwise.
	*/
	edgematch1 = (int *) ckalloc(num_endpath * sizeof(int));
	
	for (i = 0; i < num_endpath;i++) num_match1[i] = 0;
	for (i = 0; i < num_endpath;i++) edgematch1[i] = 0;

	for(i = 0; i < vertex -> num_nextedge; i ++)	{
		edge = vertex -> nextedge[i];
		/*
		  if x and y' are not connected, num_midforpath <- 0
		*/
		num_midforpath = collect2forpaths(vertex, edge1, edge, midforpath, path, num_path);  /* Px,y' */
		/*
			printf("midforpath from %d to %d\n", edge1->index, edge->index);
			PrintPaths(midforpath, num_midforpath);
		*/
		if(num_midforpath == 0)	{
		  /* 
		     If there are no paths that end (edge1,vertex->nextedge[i]), there is
		     no separating that needs to be done, so continue.
		  */
			continue;
		}
		/*	
			check consistency of P->x and Px,y'	q is not used?
		*/
		for(j = 0; j < num_endpath; j ++)	{
			c = chk_consist(&endpath[j], midforpath, num_midforpath, &q);
			/*	
			  contained <--> c = 0;	consistent <--> c-1 = # of consistent;	inconsistent <--> c = -1	
			*/
			if(c >= 0)	{
			  /*
			    Px,y' and P->x are consistent	
			  */
				if (c == 0 && (num_cont[j] == 0 || edge == edge2))	{
				  /* endpath 'j', that ends at edge1, is consistent a midforpath, a path that contains 'edge1,edge2'*/
					num_cont[j] = 1;
					if(edge == edge2)	{
					  /* The path pertains to the edge we would like to detach.*/
						edgematch1[j] = 1;
					} else	{
					  /* The path exits a different edge. */
						edgematch1[j] = 0;
					}
					num_match1[j] = 1;
				} else if(c > 0 && num_cont[j] == 0)	{
				  /* The endpath 'j' contains midforpaths that have 'edge1,edge2'.
				     This is ok only if there is just one edge that leaves 'vertex'
				     that is in the midforpath set.  If there are multiple edges leaving 'vertex'
				     then num_cont[j] will be 1, and so edgematch1[j] will remain 0.
				  */
					num_match1[j] ++;
					if(edge == edge2)	{
						edgematch1[j] = 1;
					}
				}
			}
		}
	}

/*	n1--number of P->x matching more than one Px,y'	*/

	n1 = 0;
	for(j = 0; j < num_endpath; j ++)	{
	  if(endpath[j].len_path > 0 &&  /* This path has to have edges for us to care about it.*/
	     num_startmatch > 1 &&       /* There have to be paths eminating from edge1. If nothing starts in edg1
																			then none of the paths that end in edge1 pose a problem.*/
	     (num_match1[j] == 0 ||      /* A path was found to be consistent with an edge *other* than 'edge2' 
																			that ends in edge1.*/
	      (num_match1[j] > 1 && edgematch1[j] == 1)) /* */
	     )	{
	    n1 ++;
	  }
	  /* 
	     If there is just one that starts in 'edge1', than all paths that end in 'edge1'
	     are consistent with it.
	  */
	  if(num_startmatch == 1)	{
	    edgematch1[j] = 1;
	  }
	}
	/* num_cont is the number of times a path is contained in other paths???*/
	for(j = 0; j < num_startpath; j ++)  num_cont[j] = 0;


/*	num_match2--number of branches consistent with Py->	*/

	num_match2 = (int *) ckalloc(num_startpath * sizeof(int));
	edgematch2 = (int *) ckalloc(num_startpath * sizeof(int));

	for (i = 0; i < num_startpath; i++) num_match2[i] = 0;
	for (i = 0; i < num_startpath; i++) edgematch2[i]  = 0;

	for(i = 0; i < vertex -> num_lastedge; i ++)	{
		edge = vertex -> lastedge[i];
/*	if x' and y are not connected, num_midaftpath <- 0		*/
		num_midaftpath = collect2aftpaths(vertex, edge, edge2, midaftpath, path, num_path);  /* Px',y */
		if(num_midaftpath == 0)	{
			continue;
		}
		l = 1;
		for(j = 0; j < num_startpath; j ++)	{
/*	check consistency of Py-> and Px',y	*/
		  /* Look to see if the path that starts at vertex and goes through edge2
		     is consistent with paths going through the vertex.
		  */
			c = chk_consist(&startpath[j], midaftpath, num_midaftpath, &q);
/*	contained <--> c = 0;	consistent <--> c-1 = # of consistent;	inconsistent <--> c = -1	*/
/*	x',y and y-> are consistent	*/
			if(c >= 0)	{
			  /* startpath[j] is contained in midaftpath */
			  /* All the same logic as before, except concerning the paths that go through edge2*/
				if(c == 0 && (num_cont[j] == 0 || edge == edge1))	{
					num_cont[j] = 1;
					if(edge == edge1)	{
						edgematch2[j] = 1;
					} else	{
						edgematch2[j] = 0;
					}
					num_match2[j] = 1;
				} else if(c > 0 && num_cont[j] == 0)	{
					num_match2[j] ++;
					if(edge == edge1)	{
						edgematch2[j] = 1;
					}
				}
			}
		}
	}
	/*	n2--number of Py-> matching more than one Px',y	*/

	n2 = 0;
	for(j = 0; j < num_startpath; j ++)	{
		if(startpath[j].len_path > 0 && num_endmatch > 1 && (num_match2[j] == 0 || (num_match2[j] > 1 && edgematch2[j] == 1)))	{
			n2 ++;
		}
		if(num_endmatch == 1)	{
			edgematch2[j] = 1;
		}
	}

	beginlist = (int *) ckalloc(2 * num_path * sizeof(int));
	endlist = (int *) ckalloc(2 * num_path * sizeof(int));
	nbegin = nend = 0;
	free((void *) num_cont);

	newedge = NULL;
	/*	printf("n1: %d n2: %d\n", n1, n2);*/
	if(n1 == 0 && n2 == 0)	{
	  /* 
	     There are no inconsistent paths. That means we can detach edge1 and edge2 from
	     the graph
	     edge1, edgeM->vertex->edge2, edgeN
	     becomes
	     edge'  (contains the path edge1->vertex->edge2)
	     and 
	     edgeM->vertex->edgeN
	  */
		
		begin = edge1 -> begin;
		end = edge2 -> end;
		derivelist(path, num_path, vertex, edge1, edge2, edgematch1, edgematch2,
			   beginlist, endlist, &nbegin, &nend);
		newedge = new_edge(vertex, begin, end, edge1, edge2, beginlist, endlist, nbegin, nend);
		reducepath(path, num_path, vertex, edge1, edge2, newedge, edgematch1, edgematch2);
		/* If all the paths that pass through vertex that go out through edge2 start
			 in one edge, the new edge replaces that, so we can delete it.  Otherwise, while 
			 the detach operation was valid, there still are paths that use edge1->vertex->edge2.
			 This ends up adding a *new* edge to the graph.
		*/
			 
		if(num_startmatch == 1 && vertex -> num_lastedge > 1)	{
			n = searcherase(begin -> nextedge, edge1, begin -> num_nextedge);
			erasenext(begin, n);
			n = searcherase(vertex -> lastedge, edge1, vertex -> num_lastedge);
			eraselast(vertex, n);
			/*free((void *) edge1);*/
			edge1->deleted = 1;
			/*			printf("DELETING (edge1): %d\n", edge1);*/
			f_edge[0] = 1;
			
		}
		/* same thing, but for paths that start in edge1.*/
		if(num_endmatch == 1 && edge2 != edge1 && vertex -> num_nextedge > 1)	{
			n = searcherase(vertex -> nextedge, edge2, vertex -> num_nextedge);
			erasenext(vertex, n);
			n = searcherase(end -> lastedge, edge2, end -> num_lastedge);
			eraselast(end, n);
			/*free((void *) edge2);*/
			edge2->deleted = 1;
			/*			printf("DELETING (edge2): %d\n", edge2);*/
			f_edge[1] = 1;
		}
	}

	free((void *) beginlist);
	free((void *) endlist);
	free((void *) num_match1);
	free((void *) edgematch1);
	free((void *) num_match2);
	free((void *) edgematch2);
	for(i = 0; i < num_endpath; i ++)	{
		free((void **) endpath[i].edge);
	}
	free((void *) endpath);
	for(i = 0; i < num_startpath; i ++)	{
		free((void **) startpath[i].edge);
	}
	free((void *) startpath);
	for(i = 0; i < num_midforpath; i ++)	{
		free((void **) midforpath[i].edge);
	}
	for(i = 0; i < num_midaftpath; i ++)	{
		free((void **) midaftpath[i].edge);
	}
	free((void *) midforpath);
	free((void *) midaftpath);

	if (newedge != NULL)
		assert(newedge->deleted == 0);
	return(newedge);
}
