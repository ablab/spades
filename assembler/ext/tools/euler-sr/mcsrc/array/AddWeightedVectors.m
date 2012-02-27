function vect=AddWeightedVectors(v1, w1, v2, w2)
% ADDWEIGHTEDVECTORS Add the two vectors v1 and v2 weighted according to
% w1, v2;
  
len = size(v1,2);

vect = zeros(1,len);
for i = 1:len
  v1nan = isnan(v1(i));
  v2nan = isnan(v2(i));
  if (v1nan == 0 && v2nan == 0) 
    % both values are defined, add their weighted values
    vect(i) = (v1(i)*w1 + v2(i)*w2)/(w1 + w2); 
  else
    % take the one weighte 
    if (v2nan == 1)
      vect(i) = v1(i);
    else 
      if (v1nan == 1)
	vect(i) = v2(i);
      else
	vect(i) = NaN;
      end
    end
  end
end

      

      
    

