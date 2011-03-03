ToDo:  Similarity with SNP exclusion
	   Store graph in memory
	   Dumping Sequence to disk
	   Coverage and delta integration.
TourBus:
	   "Paired" tourbus without double stranding vs."Paired tourbus" with double stranding vs. classic deBruijn(de Bruijn project) tourbus, forgetting paired structure.
		
ToThink:
	   Usage of distance: on stage of graph generation our recalculated dist? 
	   Recalculated dists after graph construction(splicing vertix using information about distance on it's neighbours) 
	   mate-reads threading?
	   
ToRefactor: 
	In: Map I from k-mers to vector of sequences.
	Out: Repeat graph for this pairs.
	
	1. Initialize Map M from (k-1) to vector (l-1)-mers. It would be our vertices. Empty for now.
	2. For each k-mers K do
	3.   For each sequence S from K.seq do
	4.     Unipath(K,S)->V_start,V_finish;
		//What is a neighbor? k-1 and k-1 intersecting on k-2 nucleo and lower sequences intersects somehow?
	5.     If V_start not present in M put it into M. 
	6.     If V_finish not present in M put it into M.
	7.     Create edge (V_start,V_finish). Store it wisely. (id of edge in vertices, sequence in file)
	8.   end do
	9.   remove (K,*) from I;
	10.end do
	
	
	Unipath(K,S):
	  K_down = K,S_down=S; 
	  while ExpandDown(K_down,S_down) is unique ExpandDown(K_down,S_down)->K_down,S_down;
	  K_up = K,S_up=S; 
	  while ExpandUp(K_up,S_up) is unique ExpandUp(K_up,S_up)->K_up,S_up;
	  V_start = (K_up,S_up);
	  V_finish = (K_down,S_down);
	return V_start, V_finish, sequence for Edge (V_start, V_finish);
	
	
	ExpandDown(K,S) 
	  Look up for possible extention to right side of (K,S) in M.
	  if find (K',S') return (not unique), (K',S')
	  else
	    look up for possible extention to right side of (K,S) in I.
		if find more than one or none return (not unique) (K,S);
		if find only one (K',S')
		  look up for possible extention to left side of (K',S')
		  if find only (K,S) return (unique), (K',S')
		  else return (not unique), (K',S')
	
	ExpandUp(K,S) similar to ExpandDown with changed direction. 
	  
	    
ToThink2: 
	(may be it is a good idea to..)
	
	1) Another ways of generating sequence from k-l mers. Now we expand k-l mer while it can be done uniquely to left, and it expansion can be uniquely expand to right and vice verse. Instead of sequence graph(?) suffix tree(?)
	2) Clustering not only from k to l (downwards), but also from sequences to kmers. (text about this idea is halfly written).
	