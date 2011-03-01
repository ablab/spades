Here is still nothing valuable!!!


ToDo:
In: Map I from k-mers to vector of sequences.
Out: Repeat graph for this pairs.

1. Initialize Map M from (k-1) to vector (l-1)-mers. It would be our vertices. Empty for now.
2. For each k-mers K do
3.   For each sequence S from K.seq do
4.     Unipath(K,S)->V_start,V_finish;
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
	  
	