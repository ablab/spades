function dist = VectorDistanceWNaN(v1, v2) 
% VECTORDISTANCEWNAN Compute the euclidian distance between two vectors
% that might have NaN values in them, interpolating the distance.
% 
%   dist = VECTORDISTANCEWNAN(v1, v2) - Given two vectors,v1 and v2, find
%   the euclidian distance between them.  Replace distances of NaN
%   components with the average component distance.  So the total
%   distance is : measured dist  + (measured dist / #non NaN) * # NaN
%   
%   v1 and v2 should be the same length
  numFinite = 0;
  numNan = 0;
  measuredDist = 0;
  ncols = size(v1,2);
  isnan1 = isnan(v1);
  isnan2 = isnan(v2);
  for i = 1:ncols
    i1p = v1(i);
    i2p = v2(i);
    if (isnan1(i) == 0 && isnan2(i) == 0)
      numFinite = numFinite + 1;
      measuredDist = measuredDist + abs(i1p - i2p);
    else
      numNan = numNan + 1;
    end
  end
  if (numNan == ncols)
    dist = NaN;
  else
    dist = measuredDist + (numNan * measuredDist / numFinite);
  end
  
  

  
  
