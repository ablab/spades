function ex =  SelectExperiment(data, cloneind, expind, expdef)
  % data is the clone-array data
  % cloneind is the index of the clone  to look for
  % expind is the experiment to consider
  % expdef is the definition of experiments, in the format
  % [start1 end1; start2 end2; start3 end3]
    
  colstart = expdef(expind,1);
  colend   = expdef(expind,2);
  orig     = data(cloneind, colstart:colend);
  ex = CleanNaN(orig);
  
