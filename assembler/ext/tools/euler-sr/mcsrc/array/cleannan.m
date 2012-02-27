function v = CleanNaN(vect)
  fi = isfinite(vect);
  nfinite = sum(fi);
  if (nfinite == 0) 
    v = [];
  else
    v = vect(fi);
  end

     
