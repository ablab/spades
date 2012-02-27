function i=FindClosestIndices(points,vals)
size(vals)
for v=1:size(vals,1)
	minP = 1;
  minDiff = abs(points(minP) - vals(v));
	for p =2:size(points,1)
		curDiff = abs(points(p) - vals(v));
		if (curDiff < minDiff)
			minDiff = curDiff;
			minP = p;
		end
	end
	i(v) = minP;
end
